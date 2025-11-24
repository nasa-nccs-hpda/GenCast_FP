import logging
import argparse
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path


def _open_xr_cf_safe(path: Path) -> xr.Dataset:
    """
    Open a NetCDF that may have illegal CF attrs (e.g., 'dtype' on time).
    We open with *no* CF decoding, strip offenders, then decode CF.
    """
    ds = xr.open_dataset(
        path,
        engine="netcdf4",
        decode_cf=False,        # <— important
        decode_times=False,     # <— important
        mask_and_scale=False,   # <— important
    )

    # Strip illegal keys from both attrs and encoding
    for name in list(ds.variables):
        if "dtype" in ds[name].attrs:
            ds[name].attrs.pop("dtype", None)
        if "dtype" in ds[name].encoding:
            ds[name].encoding.pop("dtype", None)
        # For CF coords, avoid fill values in encoding
        if name == "time":
            ds[name].encoding.pop("_FillValue", None)

    # Now decode CF safely
    ds = xr.decode_cf(ds, use_cftime=False)
    return ds


def _resolve_dt(ctime, ref_date):
    """
    Return a pandas.Timestamp for the step time.
    - If ctime is absolute (datetime64), use it directly.
    - If ctime is a timedelta64, add it to ref_date.
    - If ctime is numeric, interpret as hours since ref_date.
    """
    c = np.asarray(ctime)
    # normalize ref_date to ns precision
    ref_ns = np.datetime64(pd.to_datetime(ref_date).to_datetime64(), 'ns')

    if np.issubdtype(c.dtype, np.datetime64):
        # already absolute
        return pd.to_datetime(c.astype('datetime64[ns]'))

    if np.issubdtype(c.dtype, np.timedelta64):
        # relative offset
        return pd.to_datetime(ref_ns + c.astype('timedelta64[ns]'))

    # fallback: numeric hours since ref_date
    return pd.to_datetime(ref_ns) + pd.to_timedelta(float(c), unit='h')


def proc_time_step(
            ds_org, ctime, ref_date,
            output_dir: Path, case="init", ens_mean=True
        ):
    
    GRAV = 9.80665
    FILL_VALUE = np.float32(1.0e15)

    ds = ds_org.sel(time=ctime).expand_dims("time")

    # Time
    # dt = pd.to_datetime(ref_date + ctime)
    # instead of: pd.to_datetime(ref_date + ctime)
    dt = _resolve_dt(ctime, ref_date)
    HH = dt.strftime("%H")
    YYYY = dt.strftime("%Y")
    MM = dt.strftime("%m")
    DD = dt.strftime("%d")
    tstamp = dt.strftime("%Y-%m-%dT%H")
    begin_date = np.int32(f"{YYYY}{MM}{DD}")
    begin_time = np.int32(dt.hour * 10000)
    time_increment = np.int32(120000)
    units = f"hours since {YYYY}-{MM}-{DD} {HH}:00:00"

    ds["time"] = np.float32((ds["time"] - ds["time"]) / np.timedelta64(1, "h"))
    ds.time.attrs = {
        "long_name": "time",
        "units": units,
        "calendar": "proleptic_gregorian",
        "begin_date": begin_date,
        "begin_time": begin_time,
        "time_increment": time_increment,
    }

    # --- lat, lon, lev, ensemble ---
    lats = ds["lat"].values
    if lats[0] > lats[-1]:
        ds = ds.sel(lat=slice(None, None, -1))
    ds.lat.attrs = {"long_name": "latitude", "units": "degrees_north"}

    lons = ds["lon"].values
    if min(lons) == 0:
        ds["lon"] = ((ds["lon"] + 180) % 360) - 180
        ds = ds.sortby(ds.lon)
    ds.lon.attrs = {"long_name": "longitude", "units": "degrees_east"}

    ds = ds.rename({"level": "lev"})
    levs = ds["lev"].values.astype(np.float32)
    ds["lev"] = levs
    if levs[0] < levs[-1]:
        ds = ds.sel(lev=slice(None, None, -1))
    ds.lev.attrs = {"long_name": "pressure_level", "units": "hPa"}

    if "sample" in ds.dims:
        ds = ds.rename({"sample": "ens"})
        ds.ens.attrs = {"long_name": "ensemble_member", "units": " "}
        if ens_mean:
            ds = ds.mean(dim="ens")

    # --- variable renames + attrs (only if present) ---
    rename_dict = {
        "10m_u_component_of_wind": "U10M",
        "10m_v_component_of_wind": "V10M",
        "2m_temperature": "T2M",
        "geopotential": "H",
        "mean_sea_level_pressure": "SLP",
        "sea_surface_temperature": "SST",
        "specific_humidity": "QV",
        "temperature": "T",
        "total_precipitation_12hr": "PRECTOT",
        "u_component_of_wind": "U",
        "v_component_of_wind": "V",
        "vertical_velocity": "OMEGA",
        "geopotential_at_surface": "PHIS",
    }
    varMap = {
        "U10M": {"long_name": "10-meter_eastward_wind", "units": "m s-1"},
        "V10M": {"long_name": "10-meter_northward_wind", "units": "m s-1"},
        "T2M":  {"long_name": "2-meter_air_temperature", "units": "K"},
        "H":    {"long_name": "height", "units": "m"},
        "SLP":  {"long_name": "sea_level_pressure", "units": "Pa"},
        "SST":  {"long_name": "sea_surface_temperature", "units": "K"},
        "QV":   {"long_name": "specific_humidity", "units": "kg kg-1"},
        "T":    {"long_name": "air_temperature", "units": "K"},
        "PRECTOT": {"long_name": "total_precipitation", "units": "m"},
        "U":    {"long_name": "eastward_wind", "units": "m s-1"},
        "V":    {"long_name": "northward_wind", "units": "m s-1"},
        "OMEGA": {
            "long_name": "vertical_pressure_velocity", "units": "Pa s-1"},
        "PHIS": {
            "long_name": "surface_geopotential_height", "units": "m+2 s-2"},
    }

    valid_rename = {k: v for k, v in rename_dict.items() if k in ds.variables}
    if valid_rename:
        ds = ds.rename(valid_rename)

    for v in ds.data_vars:
        if v in varMap:
            ds[v].attrs = {
                **varMap[v],
                "_FillValue": FILL_VALUE,
                "missing_value": FILL_VALUE,
                "fmissing_value": FILL_VALUE,
            }
    
    # --- geopotential to height ---
    ds['H'] = ds['H']/GRAV

    # --- globals ---
    ds.attrs = {
        # TODO: FIX THIS TIME with +12?
        "title": f"FMGenCast forecast start at {YYYY}-{MM}-{DD}T{HH}:00:00",
        "institution": "NASA CISTO Data Science Group",
        "source": "FMGenCast model output",
        "Conventions": "CF",
        "Comment": "NetCDF-4",
    }

    # --- write ---
    compression = {"zlib": True, "complevel": 1, "shuffle": True}
    encoding = {var: compression for var in ds.data_vars}

    if case == "init":
        fname = f"FMGenCast-initial-geos_date-{tstamp}_res-1.0_levels-13.nc"
    else:
        suffix = "_ens-mean.nc" if ens_mean else ".nc"
        fname = \
            "FMGenCast-prediction-geos_date-" + \
            f"{tstamp}_res-1.0_levels-13{suffix}"

    output_dir.mkdir(parents=True, exist_ok=True)

    # removing batch dimension during postprocessing
    ds.isel(batch=0).to_netcdf(
        output_dir / fname,
        encoding=encoding,
        engine="netcdf4"
    )
    return


def run_postprocess_day(
                geos_dir: str,
                pred_dir: str,
                post_out_dir: str,
                date: str,
                ens_mean: bool = True
            ) -> None:
    """
    Process one day's init (from GEOS) and
    prediction files into CF NetCDFs.
    """
    geos_dir = Path(geos_dir)
    pred_dir = Path(pred_dir)

    Y = date.year
    M = date.month
    D = date.day
    H = date.hour

    # setup output directory
    out_day = \
        Path(post_out_dir) / f"Y{Y:04d}" / \
        f"M{M:02d}" / f"D{D:02d}" / f"T{H:02d}"

    # make output directory
    out_day.mkdir(parents=True, exist_ok=True)

    # Initial conditions (first two steps)
    init_files = sorted(
        geos_dir.glob(f"*source-geos*{Y:04d}-{M:02d}-{D:02d}T{H:02d}_*.nc"))

    if init_files:
        # ds_init = xr.open_dataset(
        # init_files[0]).drop_vars("land_sea_mask", errors="ignore")
        ds_init = _open_xr_cf_safe(init_files[0]).drop_vars(
            "land_sea_mask", errors="ignore")
        # ref_init = np.datetime64(f"{Y}-{M}-{D}T00:00:00")
        ref_init = date  # pd.Timestamp(f"{Y}-{M}-{D}T{H:02d}:00:00")
        for ctime in ds_init.time.values[:2]:
            proc_time_step(
                ds_init, ctime, ref_init,
                output_dir=out_day,
                case="init", ens_mean=ens_mean
            )
    else:
        logging.warning(
            f"No GEOS init files found for {Y}-{M}-{D}:{H:02d} in {geos_dir}")

    # Predictions (all steps)
    pred_files = sorted(
        pred_dir.glob(f"*geos_date-{Y:04d}-{M:02d}-{D:02d}T{H:02d}_*.nc"))
    if pred_files:
        # ds_pred = xr.open_dataset(
        # pred_files[0]).drop_vars("land_sea_mask", errors="ignore")
        ds_pred = _open_xr_cf_safe(
            pred_files[0]).drop_vars("land_sea_mask", errors="ignore")
        # ref_pred = np.datetime64(f"{Y}-{M}-{D}T12:00:00")
        # pd.Timestamp(f"{Y}-{M}-{D}T12:00:00") # TODO: Modify to be +12?
        ref_pred = date
        for ctime in ds_pred.time.values:
            proc_time_step(
                ds_pred, ctime, ref_pred,
                output_dir=out_day, case="pred",
                ens_mean=ens_mean
            )
    else:
        logging.warning(
            f"No prediction files found for {Y}-{M}-{D}:{H:02d} in {pred_dir}")

    return


def run_postprocess_multiday(
    start_date: str,
    end_date: str,
    geos_dir: str,
    pred_dir: str,
    post_out_dir: str,
    ens_mean: bool = True,
):
    """
    Postprocess multiple days (inclusive) of
    GenCast outputs into CF-compliant NetCDFs.
    Calls run_postprocess_day for each day in [start_date, end_date].
    """
    # start_date = np.datetime64(start_date)
    # end_date   = np.datetime64(end_date)
    # date_range = np.arange(
    # start_date, end_date + np.timedelta64(1, "D"),
    # dtype="datetime64[D]")
    fmt = "%Y-%m-%d:%H"

    # Parse exact hour from input
    start_ts = pd.to_datetime(start_date, format=fmt)
    end_ts = pd.to_datetime(end_date,   format=fmt)

    # Generate a date range in 12-hour increments
    date_range = pd.date_range(start=start_ts, end=end_ts, freq="12h")

    for current_date in date_range:

        logging.info("======================================================")
        logging.info(f"Postprocessing date: {current_date}")
        run_postprocess_day(
            geos_dir=geos_dir,
            pred_dir=pred_dir,
            post_out_dir=post_out_dir,
            date=current_date,
            ens_mean=ens_mean,
        )
        logging.info("Done postprocessing.")
        logging.info("======================================================")

    return post_out_dir


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(
        description="Convert GenCast outputs to CF-compliant NetCDFs")
    parser.add_argument(
        "--start_date", type=str, required=True,
        help="Start date (YYYY-MM-DD:HH)")
    parser.add_argument(
        "--end_date",   type=str, required=True,
        help="End date (YYYY-MM-DD:HH)")
    parser.add_argument(
        "--geos_dir", type=str, required=True,
        help="Directory with GEOS inputs (for initial conditions)")
    parser.add_argument("--pred_dir",   type=str, required=True,
                        help="Directory with GenCast predictions")
    parser.add_argument(
        "--post_out_dir", type=str, default="./output/postprocess",
        help="Directory for CF-compliant NetCDF outputs")
    parser.add_argument(
        "--no_ens_mean", action="store_true",
        help="Disable ensemble mean (keep all ensemble members)")

    args = parser.parse_args()

    run_postprocess_multiday(
        start_date=args.start_date,
        end_date=args.end_date,
        geos_dir=args.geos_dir,
        pred_dir=args.pred_dir,
        post_out_dir=args.post_out_dir,
        ens_mean=not args.no_ens_mean,
    )
