import logging
import argparse
import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path


def proc_time_step(
    ds_org, ctime, ref_date, output_dir: Path, case="init", ens_mean=True
):
    FILL_VALUE = np.float32(1.0e15)
    fmodel = "FMGenCast"

    ds = ds_org.sel(time=ctime).expand_dims("time")

    # --- time attrs ---
    # Convert to pandas timestamp for easier datetime operations
    dt = pd.Timestamp(ctime)

    # Format datetime components
    HH = f"{dt.hour:02d}"
    YYYY = f"{dt.year:04d}"
    MM = f"{dt.month:02d}"
    DD = f"{dt.day:02d}"
    tstamp = f"{YYYY}-{MM}-{DD}T{HH}"

    # Create time attributes
    begin_date = np.int32(f"{YYYY}{MM}{DD}")
    begin_time = np.int32(dt.hour * 10000)
    time_increment = np.int32(120000)
    units = f"hours since {YYYY}-{MM}-{DD} {HH}:00:00"

    # Set time coordinate to 0.0 (hours since reference time)
    # This makes the time coordinate xarray-compatible
    ds = ds.assign_coords(time=[np.float32(0.0)])
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
        "T2M": {"long_name": "2-meter_air_temperature", "units": "K"},
        "H": {"long_name": "height", "units": "m"},
        "SLP": {"long_name": "sea_level_pressure", "units": "Pa"},
        "SST": {"long_name": "sea_surface_temperature", "units": "K"},
        "QV": {"long_name": "specific_humidity", "units": "kg kg-1"},
        "T": {"long_name": "air_temperature", "units": "K"},
        "PRECTOT": {"long_name": "total_precipitation", "units": "m"},
        "U": {"long_name": "eastward_wind", "units": "m s-1"},
        "V": {"long_name": "northward_wind", "units": "m s-1"},
        "OMEGA": {
            "long_name": "vertical_pressure_velocity",
            "units": "Pa s-1",
        },
        "PHIS": {
            "long_name": "surface_geopotential_height",
            "units": "m+2 s-2",
        },
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

    # --- globals ---
    ds.attrs = {
        "title": f"FMGenCast forecast start at {YYYY}-{MM}-{DD}T12:00:00",
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
        fname = f"FMGenCast-prediction-geos_date-{tstamp}_res-1.0_levels-13{suffix}"

    output_dir.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(output_dir / fname, encoding=encoding, engine="netcdf4")


def process_dataset_files(
    file_pattern: str,
    search_dir: Path,
    ref_time_str: str,
    case_name: str,
    output_dir: Path,
    ens_mean: bool,
    Y: str,
    M: str,
    D: str,
    max_steps: int = None,
) -> None:
    """
    Process dataset files for a given pattern and case (init, pred).
    """
    files = sorted(search_dir.glob(file_pattern))
    if not files:
        logging.warning(
            f"No {case_name} files found for {Y}-{M}-{D} in {search_dir}"
        )
        return

    # Safely open dataset with time attribute handling
    try:
        # Try normal opening first
        ds = xr.open_dataset(files[0]).drop_vars(
            "land_sea_mask", errors="ignore"
        )
    except ValueError as e:
        if "dtype in attrs on variable 'time'" in str(e):
            logging.info("Handling time dtype attribute conflict...")
            # Open without time decoding first
            ds = xr.open_dataset(files[0], decode_times=False).drop_vars(
                "land_sea_mask", errors="ignore"
            )

            # Clean problematic time attributes
            if "time" in ds.variables:
                # Store original attributes we might need
                time_attrs = dict(ds.time.attrs)

                # Remove problematic encoding attributes
                problematic_attrs = ["dtype", "_FillValue", "missing_value"]
                for attr in problematic_attrs:
                    if attr in ds.time.attrs:
                        ds.time.attrs.pop(attr)

                # Try to decode times manually
                try:
                    ds = xr.decode_cf(ds, decode_times=True)
                    logging.info(
                        "Successfully decoded times after cleaning attributes"
                    )
                except Exception as decode_error:
                    logging.warning(f"Could not decode times: {decode_error}")
                    # Restore some attributes if decoding failed
                    ds.time.attrs.update(
                        {
                            k: v
                            for k, v in time_attrs.items()
                            if k not in problematic_attrs
                        }
                    )
        else:
            raise

    # Create reference time as datetime64
    ref_time = np.datetime64(f"{Y}-{M}-{D}T{ref_time_str}")

    # Process time steps - handle both decoded and non-decoded time cases
    time_values = ds.time.values[:max_steps] if max_steps else ds.time.values

    for ctime in time_values:
        try:
            # Ensure ctime is datetime64
            if isinstance(ctime, (int, float, np.integer, np.floating)):
                # If time is numeric, try to convert using time attributes
                if "units" in ds.time.attrs:
                    # Use pandas/xarray to convert numeric time to datetime
                    ctime_dt64 = pd.to_datetime(
                        ctime,
                        unit="D",
                        origin=ds.time.attrs.get(
                            "units", "days since 1900-01-01"
                        ),
                    )
                    ctime_dt64 = np.datetime64(ctime_dt64)
                else:
                    logging.warning(
                        f"Numeric time value {ctime} without units - skipping"
                    )
                    continue
            else:
                # Already a datetime-like object
                ctime_dt64 = np.datetime64(ctime)

            proc_time_step(
                ds,
                ctime_dt64,
                ref_time,
                output_dir=output_dir,
                case=case_name,
                ens_mean=ens_mean,
            )
        except Exception as time_error:
            logging.error(f"Error processing time step {ctime}: {time_error}")
            continue


def run_postprocess_day(
    geos_dir: str,
    pred_dir: str,
    post_out_dir: str,
    year: int,
    month: int,
    day: int,
    ens_mean: bool = True,
) -> None:
    """Process one day's init (from GEOS) and prediction files into CF NetCDFs."""
    geos_dir = Path(geos_dir)
    pred_dir = Path(pred_dir)
    out_day = (
        Path(post_out_dir) / f"Y{year:04d}" / f"M{month:02d}" / f"D{day:02d}"
    )
    out_day.mkdir(parents=True, exist_ok=True)

    # Format date strings
    Y = f"{year:04d}"
    M = f"{month:02d}"
    D = f"{day:02d}"

    # Hard-coded configuration for processing
    INIT_CONFIG = {
        "file_pattern": f"*source-geos*{Y}-{M}-{D}_*.nc",
        "ref_time_str": "00:00:00",
        "case_name": "init",
        "max_steps": 2,  # First two steps only
    }

    PRED_CONFIG = {
        "file_pattern": f"*geos_date-{Y}-{M}-{D}_*.nc",
        "ref_time_str": "12:00:00",
        "case_name": "pred",
        "max_steps": None,  # All steps
    }

    # Process initial conditions
    process_dataset_files(
        INIT_CONFIG["file_pattern"],
        geos_dir,
        INIT_CONFIG["ref_time_str"],
        INIT_CONFIG["case_name"],
        out_day,
        ens_mean,
        Y,
        M,
        D,
        INIT_CONFIG["max_steps"],
    )

    # Process predictions
    process_dataset_files(
        PRED_CONFIG["file_pattern"],
        pred_dir,
        PRED_CONFIG["ref_time_str"],
        PRED_CONFIG["case_name"],
        out_day,
        ens_mean,
        Y,
        M,
        D,
        PRED_CONFIG["max_steps"],
    )


def run_postprocess_multiday(
    start_date: str,
    end_date: str,
    geos_dir: str,
    pred_dir: str,
    post_out_dir: str,
    ens_mean: bool = True,
):
    """Postprocess multiple days (inclusive) of GenCast outputs into CF-compliant NetCDFs.
    Calls run_postprocess_day for each day in [start_date, end_date].
    """
    start_date = np.datetime64(start_date)
    end_date = np.datetime64(end_date)
    date_range = np.arange(
        start_date, end_date + np.timedelta64(1, "D"), dtype="datetime64[D]"
    )

    for current_date in date_range:
        y = int(str(current_date)[:4])
        m = int(str(current_date)[5:7])
        d = int(str(current_date)[8:10])

        logging.info("======================================================")
        logging.info(f"Postprocessing date: {current_date}")
        run_postprocess_day(
            geos_dir=geos_dir,
            pred_dir=pred_dir,
            post_out_dir=post_out_dir,
            year=y,
            month=m,
            day=d,
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
        description="Convert GenCast outputs to CF-compliant NetCDFs"
    )
    parser.add_argument(
        "--start_date", type=str, required=True, help="Start date (YYYY-MM-DD)"
    )
    parser.add_argument(
        "--end_date", type=str, required=True, help="End date (YYYY-MM-DD)"
    )
    parser.add_argument(
        "--geos_dir",
        type=str,
        required=True,
        help="Directory with GEOS inputs (for initial conditions)",
    )
    parser.add_argument(
        "--pred_dir",
        type=str,
        required=True,
        help="Directory with GenCast predictions",
    )
    parser.add_argument(
        "--post_out_dir",
        type=str,
        default="./output/postprocess",
        help="Directory for CF-compliant NetCDF outputs",
    )
    parser.add_argument(
        "--no_ens_mean",
        action="store_true",
        help="Disable ensemble mean (keep all ensemble members)",
    )

    args = parser.parse_args()

    run_postprocess_multiday(
        start_date=args.start_date,
        end_date=args.end_date,
        geos_dir=args.geos_dir,
        pred_dir=args.pred_dir,
        post_out_dir=args.post_out_dir,
        ens_mean=not args.no_ens_mean,
    )
