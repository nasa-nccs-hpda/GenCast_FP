# This is the driver code for converting GEOS-FP data to ERA5 format

import os
import argparse
import logging
import numpy as np
import xesmf as xe
import xarray as xr
import pandas as pd

# from fp_to_era5 import *
from gencast_fp.preprocess.fp_to_era5 import (
    discover_files,
    fp_to_era5_xlevs,
    era5_dataset,
    sst_dataset,
    fp_to_era5_hgrid,
)

def generate_dates(start_date: str, end_date: str):
    fmt = "%Y-%m-%d:%H"
    try:
        start = pd.to_datetime(start_date, format=fmt)
        end = pd.to_datetime(end_date, format=fmt)
    except:
        raise ValueError("Please provide dates in format: YYYY-MM-DD-HH (e.g., '2020-01-01-00')")
    start_shift = start - pd.Timedelta(hours=12)
    dates = pd.date_range(start=start_shift, end=end, freq="12h")
    return dates

def get_era5_lsm(lsm_file: str = "/css/era5/static/era5_static-allvar.nc"):
    ds = xr.open_dataset(lsm_file, engine="netcdf4")
    lsm = (
        ds["lsm"]
        .squeeze(drop=True)
        .drop_vars(["expver", "number"])
        .astype("float32")
    )
    return lsm

def get_sst(sst_file: str, current_date: pd.Timestamp):
    # Read the binary SST file and extract the SST for the given date
    # The SST file is in a custom binary format with a header
    # The header contains the start date of the record
    # Each record is 72 bytes + 2880*1440*4 bytes (float32)
    # The first 68 bytes of the header are not used
    nx, ny = 2880, 1440
    dtype = np.float32
    theader_0 = 68
    header_count = theader_0 // 4
    theader = 72
    # Get the start_date (yyyy-mm-dd) of the SST record
    with open(sst_file, "rb") as f:
        header = np.fromfile(f, dtype=dtype, count=header_count)
    rec_start_date = pd.Timestamp(int(header[1]), int(header[2]), int(header[3]))
    delta = (current_date - rec_start_date).days
    if delta < 0:
        raise ValueError(f"SST file {sst_file} does not cover date {current_date}")
    
    # Read the SST record for the current date
    offset = theader_0 + delta * (theader + nx * ny * 4)
    with open(sst_file, "rb") as f:
        f.seek(offset, 0)
        data = np.fromfile(f, dtype=dtype, count=nx * ny)
        sst = data.reshape((ny, nx))
    
    # Convert to xarray DataArray
    ds = sst_dataset("OSTIA-REYNOLDS on ERA-5 Grid for AI/ML Modeling")
    sst_np_time = sst[np.newaxis, :, :]  # add time dimension
    time = xr.DataArray([current_date], dims=["time"], coords={"time": [current_date]})
    ds = ds.assign_coords(time=time)
    ds['sst'] = xr.DataArray( sst_np_time, 
                              dims=("time", "latitude", "longitude"), 
                              coords={"time": time, "latitude": ds['latitude'], "longitude": ds['longitude']},
                              attrs={"units": "K", 
                                     "long_name": "Sea Surface Temperature",
                                    "standard_name": "sea_surface_temperature"},
                            )

    return ds

def expand_dims(ds, steps):
    # Expand the time dimension of the dataset
    orig_time = ds.time.values
    extra_steps = steps - 2

    ds_last = ds.isel(time=1)
    new_times = pd.date_range(
        start=orig_time[-1], periods=extra_steps + 1, freq="12h"
    )[1:]

    repeated = ds_last.expand_dims(time=range(extra_steps)).copy(deep=True)
    repeated["time"] = new_times
    ds_extended = xr.concat([ds, repeated], dim="time")
    return ds_extended


def to_gencast_input(ds):

    nsteps = 32
    GRAV = 9.80665

    levs = np.array(
        [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 850, 925, 1000]
    )

    static = [
        "land_sea_mask",
        "geopotential_at_surface",
    ]

    var_2d = [
        "2m_temperature",
        "sea_surface_temperature",
        "mean_sea_level_pressure",
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "total_precipitation",
    ]

    var_3d = [
        "temperature",
        "specific_humidity",
        "u_component_of_wind",
        "v_component_of_wind",
        "vertical_velocity",
        "geopotential",
    ]
    var_mapping = {
        "t2m": "2m_temperature",
        "t": "temperature",
        "u10": "10m_u_component_of_wind",
        "v10": "10m_v_component_of_wind",
        "u": "u_component_of_wind",
        "v": "v_component_of_wind",
        "q": "specific_humidity",
        "w": "vertical_velocity",
        "z": "geopotential",
        "msl": "mean_sea_level_pressure",
        "tp": "total_precipitation_12hr",
        "zs": "geopotential_at_surface",
        "sst": "sea_surface_temperature",
        "lsm": "land_sea_mask",
        "latitude": "lat",
        "longitude": "lon",
        "pressure_level": "level",
    }

    # coarsen the data & reverse the latitude
    ds = ds.isel(
        latitude=slice(None, None, -4), longitude=slice(None, None, 4)
    ).compute()

    ds = ds.drop_vars(["hgt", "p", "sp", "skt"])
    # change variable names
    ds = ds.rename(var_mapping)

    # change data types
    ds = ds.astype({var: "float32" for var in ds.data_vars})

    # expand the time dimension
    lsm = ds["land_sea_mask"]
    ds = ds.drop_vars(["land_sea_mask"])
    ds = expand_dims(ds, nsteps)

    ds = ds.assign_coords(datetime=ds["time"])

    # expand the dimensions
    for var in list(ds.data_vars) + ["datetime"]:
        ds[var] = ds[var].expand_dims("batch")

    # TODO: Commenting it for now to see if GenCast complains
    # change time coordinate to timedelta
    # ds["time"] = ds["time"] - ds["time"].isel(time=0)

    # # add land_sea_mask
    # file = f"/discover/nobackup/projects/QEFM/data/FMGenCast/12hr/Y2024/gencast-dataset-source-era5_date-{date_str}_res-1.0_levels-13_steps-20.nc"
    # ds_lsm = xr.open_dataset(file)
    # ds['land_sea_mask'] = ds_lsm['land_sea_mask']
    # # using the sea_surface_temperature from the original dataset
    # ds['sea_surface_temperature'] = ds_lsm['sea_surface_temperature']

    # drop the time dimension for geopotential_at_surface
    ds["geopotential_at_surface"] = (
        ds["geopotential_at_surface"].isel(time=0).drop_vars(["time"])
    )
    ds["land_sea_mask"] = lsm

    return ds


def run_preprocess(start_date, end_date, outdir, expid):

    os.makedirs(outdir, exist_ok=True)

    res_value = 1.0  # 1.0 resolution
    nsteps = 30  # 15 day rollout

    dates = generate_dates(start_date, end_date)
    pairs = [(dates[i], dates[i+1]) for i in range(len(dates)-1)]

    regridder = None
    sst_regridder = None
    for time_tuple in pairs:

        date_str = time_tuple[1].strftime("%Y-%m-%dT%H")
        out_file = os.path.join(
            outdir,
            f"gencast-dataset-source-geos_date-{date_str}"
            + f"_res-{res_value}_levels-13_steps-{nsteps}.nc",
        )

        # skip if already exists
        if os.path.exists(out_file):
            logging.info(f"Skipping {out_file}, already exists.")
            continue

        daily_Ex = []
        daily_Ep = []
        for dt in time_tuple:
            logging.info(dt)
            Files = discover_files(dt, outdir=outdir, expid=expid)
            logging.info(Files)

            fp_Nx = xr.open_dataset(Files["fp_Nx"], engine="netcdf4")
            fp_Nv = xr.open_dataset(Files["fp_Nv"], engine="netcdf4")

            ai_Nx, ai_Np = fp_to_era5_xlevs(fp_Nx, fp_Nv)

            ai_Ex = era5_dataset(
                "GEOS-FP 2D Variables on ERA-5 Grid for AI/ML Modeling"
            )

            if not regridder:
                regridder = xe.Regridder(ai_Nx, ai_Ex, "conservative")

            ai_Ex, ai_Ep = fp_to_era5_hgrid(ai_Nx, ai_Np, regridder=regridder)

            # get the sst
            sst_ds = get_sst(Files["sst"], dt)
            print(sst_ds['sst'].min(), sst_ds['sst'].max())
            exit()

            if not sst_regridder:
                sst_grid = sst_dataset("OSTIA-REYNOLDS on ERA-5 Grid for AI/ML Modeling")
                sst_regridder = xe.Regridder(sst_grid, ai_Ex, "conservative")
            ai_sst = sst_regridder(sst_ds['sst'], keep_attrs=True)
            ai_Ex['sst'] = ai_sst

            daily_Ex.append(ai_Ex)
            daily_Ep.append(ai_Ep)

        # concat along time
        ai_Ex_day = xr.concat(daily_Ex, dim="time")
        ai_Ep_day = xr.concat(daily_Ep, dim="time")

        # add the lsm
        lsm = get_era5_lsm(Files["e5_Es"])

        # apply lsm to sst
        lsm_nan = lsm.where(lsm == 0, 1.0, np.nan)
        ai_Ex_day['sst'] = ai_Ex_day['sst'] * lsm_nan

        # merge into single dataset for the day
        ai_day = xr.merge([ai_Ex_day, ai_Ep_day, lsm.to_dataset()])

        ds_out = to_gencast_input(ai_day)

        # save to netcdf
        # date_str = dt.strftime("%Y-%m-%d")
        # out_file = (
        #    f"{outdir}/gencast-dataset-source-geos_date-{date_str}"
        #    f"_res-{res_value}_levels-13_steps-{nsteps}.nc"
        # )
        ds_out.to_netcdf(
            out_file, mode="w", format="NETCDF4", engine="netcdf4"
        )
        exit()
    return


def main():

    parser = argparse.ArgumentParser(
        description="Convert GEOS-FP data to ERA5 format."
    )
    parser.add_argument(
        "--start_date",
        type=str,
        required=True,
        help="Start date to process (YYYY-MM-DD:HH)",
    )
    parser.add_argument(
        "--end_date",
        type=str,
        required=True,
        help="End date to process (YYYY-MM-DD:HH)",
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default="./output/",
        help="Output directory for the converted files",
    )
    parser.add_argument(
        "--expid",
        type=str,
        default="f5295",
        help="Experiment ID for the output files",
    )

    args = parser.parse_args()

    # Run preprocessing function
    run_preprocess(args.start_date, args.end_date, args.outdir, args.expid)

    return


if __name__ == "__main__":
    main()
