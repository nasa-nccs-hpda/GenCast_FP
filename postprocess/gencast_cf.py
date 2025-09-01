import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import argparse

def proc_time_step(ds_org, ctime, ref_date, output_dir, case="init"):
    fmodel = "FMGenCast"
    FILL_VALUE = np.float32(1.e+15)
    ds = ds_org.sel(time=ctime).expand_dims("time")


    # Time
    dt = pd.to_datetime(ref_date + ctime)
    HH = dt.strftime("%H")
    YYYY = dt.strftime("%Y")
    MM = dt.strftime("%m")
    DD = dt.strftime("%d")

    tstamp = dt.strftime("%Y-%m-%dT%H")     

    long_name = "time"
    begin_date = np.int32(f"{YYYY}{MM}{DD}")
    begin_time = np.int32(dt.hour * 10000)
    time_increment = np.int32(120000)
    units = f"hours since {YYYY}-{MM}-{DD} {HH}:00:00"
    calendar = "proleptic_gregorian"

    # change value of time
    ds['time'] = np.float32((ds['time']- ds['time'])/np.timedelta64(1, 'h'))
    # add attributes
    ds.time.attrs = {
        "long_name" : long_name,
        "units" : units,
        "calendar" : calendar,
        "begin_date" : begin_date,
        "begin_time" : begin_time,
        "time_increment": time_increment
    }

    # Latitude
    lats = ds['lat'].values
    if lats[0] > lats[-1]:
        ds = ds.sel(lat=slice(None, None, -1))

    ds.lat.attrs = {
        "long_name" : "latitude",
        "units" : "degrees_north",
    }

    # Longitude
    lons = ds['lon'].values
    if min(lons) == 0:
        ds['lon'] = ((ds["lon"] + 180) % 360) - 180
        ds = ds.sortby(ds.lon)
    ds.lon.attrs = {
        "long_name" : "longitude",
        "units" : "degrees_east",
    }

    # level
    ds = ds.rename({'level': 'lev'})
    levs = ds['lev'].values.astype(np.float32)
    ds['lev'] = levs
    if levs[0] < levs[-1]:
        # Flip the level array
    #    ds['lev'] = levs[::-1]
        # Flip the data array
        ds = ds.sel(lev=slice(None, None, -1))
    ds.lev.attrs = {
        "long_name" : "pressure_level",
        "units" : "hPa",
    }

    # ensemble
    if "sample" in ds.dims:
        ds = ds.rename({'sample': 'ens'})
        ds.ens.attrs = {
            "long_name" : "ensemble_member",
            "units" : " ",
        }

        ## Calculate ensemble mean
        if ens_mean:
            ds = ds.mean(dim="ens")

    ## Variables
    # rename variables
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
    valid_rename_dict = {k: v for k, v in rename_dict.items() if k in ds.variables}
    ds = ds.rename(valid_rename_dict)

    # map attributes
    varMap = {
        "U10M": {
            "long_name" : "10-meter_eastward_wind",
            "units" : "m s-1",
        },
        "V10M": {
            "long_name" : "10-meter_northward_wind",
            "units" : "m s-1",
        },
        "T2M": {
            "long_name" : "2-meter_air_temperature",
            "units" : "K",
        },
        "H": {
            "long_name" : "height",
            "units" : "m",
        },
        "SLP": {
            "long_name" : "sea_level_pressure",
            "units" : "Pa",
        },
        "SST": {
            "long_name" : "sea_surface_temperature",
            "units" : "K",
        },
        "QV": {
            "long_name" : "specific_humidity",
            "units" : "kg kg-1",
        },
        "T": {
            "long_name" : "air_temperature",
            "units" : "K",
        },
        "PRECTOT": {
            "long_name" : "total_precipitation",
            "units" : "m",
        },
        "U": {
            "long_name" : "eastward_wind",
            "units" : "m s-1",
        },
        "V": {
            "long_name" : "northward_wind",
            "units" : "m s-1",
        },
        "OMEGA": {
            "long_name" : "vertical_pressure_velocity",
            "units" : "Pa s-1",
        },
        "PHIS": {
            "long_name" : "surface_geopotential_height",
            "units" : "m+2 s-2",
        },
    }

    for var in ds.data_vars:
    #    # add attributes
        ds[var].attrs = varMap[var]
        ds[var].attrs['_FillValue'] = FILL_VALUE
        ds[var].attrs['missing_value'] = FILL_VALUE
        ds[var].attrs['fmissing_value'] = FILL_VALUE


    ## add global attributes
    ds.attrs = {
        "title" : f"{fmodel} forecast start at {YYYY}-{MM}-{DD}T12:00:00", 
        "institution" : "NASA CISTO Data Science Group",
        "source" : f"{fmodel} model output",
        "Conventions" : "CF",
        "Comment" : "NetCDF-4" 
    }


    # ## Write to NetCDF
    compression = {"zlib": True, 
                "complevel": 1,
                "shuffle": True,}
    encoding = {var: compression for var in ds.data_vars}

    if case == "init":
        fname = f"{fmodel}-initial-geos_date-{tstamp}_res-1.0_levels-13.nc"
    else:
        fname = f"{fmodel}-prediction-geos_date-{tstamp}_res-1.0_levels-13_ens-mean.nc"
    output_file = output_dir / fname
    ds.to_netcdf(output_file, encoding=encoding, engine="netcdf4")



    print("Finished \n", ds)

def main():

    parser = argparse.ArgumentParser(description="Convert GenCast output to CF-compliant NetCDF")
    # parser.add_argument("input_dir", type=str, help="Path to GenCast output directory")
    # parser.add_argument("fmodel", type=str, help="Model name")
    parser.add_argument("--geos_dir", "-g", type=str, help="Inputs for GenCast to get surface geopotential height")
    parser.add_argument("--pred_dir", "-p", type=str, help="directory of GenCast prediction files")
    parser.add_argument("--output_dir", "-o", type=str, help="Output directory to CF-compliant NetCDF")
    parser.add_argument("--year", "-y", type=str, help="Year")
    parser.add_argument("--month", "-m", type=str, help="Month")
    parser.add_argument("--day", "-d", type=str, help="Day")

    args = parser.parse_args()

    output_dir = Path(args.output_dir)/f"Y{args.year}"/f"M{args.month}"/f"D{args.day}"
    output_dir.mkdir(parents=True, exist_ok=True)

    for case in ["init", "pred"]:
        if case == "init":
            ## Process initial conditions
            # Retrieve files from input dir for GenCast
            files = sorted(Path(args.geos_dir).glob(f"*source-geos*{args.year}-{args.month}-{args.day}_*.nc"))
            ref_date = np.datetime64(f"{args.year}-{args.month}-{args.day}T00:00:00")
            n=2 # only process first two time steps for initial conditions
        else:
            files = sorted(Path(args.pred_dir).glob(f"*geos_date-{args.year}-{args.month}-{args.day}_*.nc"))
            ref_date = np.datetime64(f"{args.year}-{args.month}-{args.day}T12:00:00")
            n = -1 # process all time steps for prediction 
    
    ds = xr.open_dataset(files[0])
    ds = ds.drop_vars("land_sea_mask", errors='ignore')
    for ctime in ds.time.values[:n]:
        print("Processing time ", f"{ctime.strftime("%Y-%m-%d")}", " for case ", case)
        proc_time_step(ds, ctime, ref_date, output_dir=output_dir, case=case)



    # for ctime in ds_init.time.values[:2]:
    #     proc_time_step(ds_init, ctime, ref_date, output_dir=output_dir, case="init")
    # exit()    

    # ref_date = np.datetime64(f"{args.year}-{args.month}-{args.day}T12:00:00")
    # #pred_dir = Path(args.pred_dir)
    # pred_dir = Path("/discover/nobackup/projects/QEFM/data/rollout_outputs/FMGenCast/raw/geos")
    # files = sorted(pred_dir.glob(f"*geos_date-{args.year}-{args.month}-{args.day}_*.nc"))
    # file = files[0]
    # ds_org = xr.open_dataset(file)
    # ## ds_org contains all time steps for the 10-day prediction

    # output_dir = Path(args.geos_dir)

    # for ctime in ds_org.time.values:
    #     proc_time_step(ds_org, ctime, ref_date, output_dir)

if __name__ == "__main__":
    ens_mean = True
    main()