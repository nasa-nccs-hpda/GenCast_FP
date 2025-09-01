import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import argparse

def proc_time_step(ds_org, ctime, ref_date, output_dir):
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
        print("After coord \n", ds)

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
    }
    ds = ds.rename(rename_dict)
    print("After rename \n ", ds)


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
    # compression = {"zlib": True, 
    #             "complevel": 1,
    #             "shuffle": True,}
    # encoding = {var: compression for var in ds.data_vars}
    # #output_dir = Path(f"/discover/nobackup/projects/QEFM/data/rollout_outputs/{fmodel}/geos-fp-interp-no-mask/Y{yyyy}/M{mm}/D{dd}")
    # #output_dir = Path(f"/discover/nobackup/projects/QEFM/data/rollout_outputs/{fmodel}/0p25/Y{yyyy}/M{mm}/D{dd}")
    # output_dir.mkdir(parents=True, exist_ok=True)
    # fname = f"GenCast-prediction-geos_date-{tstamp}_res-1.0_levels-13_ens-mean.nc"
    # output_file = output_dir / fname
    # ds.to_netcdf(output_file, encoding=encoding, engine="netcdf4")



    print("Finished \n", ds)

def main():

    parser = argparse.ArgumentParser(description="Convert GenCast output to CF-compliant NetCDF")
    # parser.add_argument("input_dir", type=str, help="Path to GenCast output directory")
    # parser.add_argument("fmodel", type=str, help="Model name")
    parser.add_argument("--geos_dir", "-g", type=str, help="Inputs for GenCast to get surface geopotential height")
    # parser.add_argument("--pred_dir", "-p", type=str, help="directory of GenCast prediction files")
    parser.add_argument("--year", "-y", type=str, help="Year")
    parser.add_argument("--month", "-m", type=str, help="Month")
    parser.add_argument("--day", "-d", type=str, help="Day")
    args = parser.parse_args()

    ## Process initial conditions
    # Retrieve files from input dir for GenCast
    geos_dir = Path(args.geos_dir)
    files = sorted(geos_dir.glob(f"*source-geos*{args.year}-{args.month}-{args.day}_*.nc"))
    file = files[0]
    ref_date = np.datetime64(f"{args.year}-{args.month}-{args.day}T00:00:00")
    ds_init = xr.open_dataset(file)
    for ctim in ds_init.time.values[:2]:
        proc_time_step(ds_init, ctime, red_date, output_dir=None)
    exit()    

    ref_date = np.datetime64(f"{args.year}-{args.month}-{args.day}T12:00:00")
    #pred_dir = Path(args.pred_dir)
    pred_dir = Path("/discover/nobackup/projects/QEFM/data/rollout_outputs/FMGenCast/raw/geos")
    files = sorted(pred_dir.glob(f"*geos_date-{args.year}-{args.month}-{args.day}_*.nc"))
    file = files[0]
    ds_org = xr.open_dataset(file)
    ## ds_org contains all time steps for the 10-day prediction

    fmodel = "FMGenCast"
    output_dir = Path(f"/discover/nobackup/projects/QEFM/data/rollout_outputs/{fmodel}/geos-fp-interp-no-mask/Y{args.year}/M{args.month}/D{args.day}")

    for ctime in ds_org.time.values:
        proc_time_step(ds_org, ctime, ref_date, output_dir)

if __name__ == "__main__":
    ens_mean = True
    main()
# ## Example usage
# #fmodel = "FMGenCast"
# #file_path = input_dir / fmodel / 'raw' / 'Y2024'
# #file_path = input_dir / fmodel / 'raw' / '0p25'
# file_path = input_dir
# yyyy = args.year
# mm = args.month
# dd = args.day

# files = sorted(file_path.glob(f"*geos_date-{yyyy}-{mm}-{dd}_*.nc"))
# file = files[0]

# MAPL_GRAV = 9.80665
# FILL_VALUE = np.float32(1.e+15)
# ds_org = xr.open_dataset(file)
# print("At Open : \n", ds_org)
# ref_date = np.datetime64(f"{yyyy}-{mm}-{dd}T12:00:00")

# ## Coordinates
# # For GenCast Only, remove "batch"
# ds_org = ds_org.squeeze(dim="batch")

# # add variable geopotential at surface
# #source = Path("/discover/nobackup/projects/QEFM/data/FMGenCast/12hr/Y2024")
# #source = Path("/discover/nobackup/projects/QEFM/data/FMGenCast/0p25")
# source = Path(args.geos_dir)
# tmp_file = list(source.glob(f"*{yyyy}-{mm}-{dd}_*.nc"))[0]
# ds_temp = xr.open_dataset(tmp_file)
# ds_org['PHIS'] = ds_temp['geopotential_at_surface']
# ens_mean = True

# for ctime in ds_org.time.values:
#     ds = ds_org.sel(time=ctime).expand_dims("time")
#     dt = pd.to_datetime(ref_date + ctime)
#     HH = dt.strftime("%H")
#     YYYY = dt.strftime("%Y")
#     MM = dt.strftime("%m")
#     DD = dt.strftime("%d")

#     tstamp = dt.strftime("%Y-%m-%dT%H") 

#     # Time
#     long_name = "time"
#     begin_date = np.int32(f"{YYYY}{MM}{DD}")
#     begin_time = np.int32(dt.hour * 10000)
#     time_increment = np.int32(120000)
#     units = f"hours since {YYYY}-{MM}-{DD} {HH}:00:00"
#     calendar = "proleptic_gregorian"


#     # change value of time
#     ds['time'] = np.float32((ds['time']- ds['time'])/np.timedelta64(1, 'h'))
#     # add attributes
#     ds.time.attrs = {
#         "long_name" : long_name,
#         "units" : units,
#         "calendar" : calendar,
#         "begin_date" : begin_date,
#         "begin_time" : begin_time,
#         "time_increment": time_increment
#     }

#     # Latitude
#     lats = ds['lat'].values
#     fill_north = False
#     fill_south = False
#     if lats[0] > lats[-1]:
#         # Flip the latitude array
#     #    ds['lat'] = lats[::-1]
#         # Flip the data array
#         ds = ds.sel(lat=slice(None, None, -1))
#     # if -90. not in ds.lat.values:
#     #     fill_south = True
#     #     lat_values = np.insert(ds.lat.values, 0, -90)
#     #     ds = ds.assign_coords(lat=lat_values)
#     # if 90. not in ds.lat.values:
#     #     fill_north = True
#     #     lat_values = np.append(ds.lat.values, 90)
#     #     ds = ds.assign_coords(lat=lat_values)

#     ds.lat.attrs = {
#         "long_name" : "latitude",
#         "units" : "degrees_north",
#     }

#     # Longitude
#     lons = ds['lon'].values
#     if min(lons) == 0:
#         ds['lon'] = ((ds["lon"] + 180) % 360) - 180
#         ds = ds.sortby(ds.lon)
#     ds.lon.attrs = {
#         "long_name" : "longitude",
#         "units" : "degrees_east",
#     }

#     # level
#     ds = ds.rename({'level': 'lev'})
#     levs = ds['lev'].values.astype(np.float32)
#     ds['lev'] = levs
#     if levs[0] < levs[-1]:
#         # Flip the level array
#     #    ds['lev'] = levs[::-1]
#         # Flip the data array
#         ds = ds.sel(lev=slice(None, None, -1))
#     ds.lev.attrs = {
#         "long_name" : "pressure_level",
#         "units" : "hPa",
#     }

#     # ensemble
#     ds = ds.rename({'sample': 'ens'})
#     ds.ens.attrs = {
#         "long_name" : "ensemble_member",
#         "units" : " ",
#     }
#     print("After coord \n", ds)

#     ## Calculate ensemble mean
#     if ens_mean:
#         ds = ds.mean(dim="ens")

#     ## Variables
#     # rename variables
#     rename_dict = {
#         "10m_u_component_of_wind": "U10M",
#         "10m_v_component_of_wind": "V10M",
#         "2m_temperature": "T2M",
#         "geopotential": "H",
#         "mean_sea_level_pressure": "SLP",
#         "sea_surface_temperature": "SST",
#         "specific_humidity": "QV",
#         "temperature": "T",
#         "total_precipitation_12hr": "PRECTOT",
#         "u_component_of_wind": "U",
#         "v_component_of_wind": "V",
#         "vertical_velocity": "OMEGA",
#     }
#     ds = ds.rename(rename_dict)
#     print("After rename \n ", ds)


#     # map attributes
#     varMap = {
#         "U10M": {
#             "long_name" : "10-meter_eastward_wind",
#             "units" : "m s-1",
#         },
#         "V10M": {
#             "long_name" : "10-meter_northward_wind",
#             "units" : "m s-1",
#         },
#         "T2M": {
#             "long_name" : "2-meter_air_temperature",
#             "units" : "K",
#         },
#         "H": {
#             "long_name" : "height",
#             "units" : "m",
#         },
#         "SLP": {
#             "long_name" : "sea_level_pressure",
#             "units" : "Pa",
#         },
#         "SST": {
#             "long_name" : "sea_surface_temperature",
#             "units" : "K",
#         },
#         "QV": {
#             "long_name" : "specific_humidity",
#             "units" : "kg kg-1",
#         },
#         "T": {
#             "long_name" : "air_temperature",
#             "units" : "K",
#         },
#         "PRECTOT": {
#             "long_name" : "total_precipitation",
#             "units" : "m",
#         },
#         "U": {
#             "long_name" : "eastward_wind",
#             "units" : "m s-1",
#         },
#         "V": {
#             "long_name" : "northward_wind",
#             "units" : "m s-1",
#         },
#         "OMEGA": {
#             "long_name" : "vertical_pressure_velocity",
#             "units" : "Pa s-1",
#         },
#         "PHIS": {
#             "long_name" : "surface_geopotential_height",
#             "units" : "m+2 s-2",
#         },
#     }

#     # define chunk size for each variable
#     # mask variables based on elevation
#     ## TODO : Need to check surface geopotential height from ERA5
#     ## Get geopotential height 
#     ds['H'] = ds['H']/MAPL_GRAV

#     # topo
#     topo = ds.PHIS.values/MAPL_GRAV
#     height = ds.H.values
#     mask = np.where(height > topo, 1, 0)
#     for var in ds.data_vars:
#     #    # add attributes
#         ds[var].attrs = varMap[var]
#         ds[var].attrs['_FillValue'] = FILL_VALUE
#         ds[var].attrs['missing_value'] = FILL_VALUE
#         ds[var].attrs['fmissing_value'] = FILL_VALUE
#         # mask 3d variables
#    #     if 'lev' in ds[var].dims:
#    #         ds[var] = ds[var].where(mask == 1, FILL_VALUE)
#     # chunk 
#     #nlats = len(ds.lat)
#     #nlons = len(ds.lon)
#     #chunks_size = {"ens": 1, "time": 1, "lev": 1, "lat": nlats, "lon": nlons}
#     #ds = ds.chunk(chunks_size)
#     print("After variable \n", ds)

#     ## add global attributes
#     ds.attrs = {
#         "title" : f"{fmodel} forecast start at {yyyy}-{mm}-{dd}T12:00:00", 
#         "institution" : "NASA CISTO Data Science Group",
#         "source" : f"{fmodel} model output",
#         "Conventions" : "CF",
#         "Comment" : "NetCDF-4" 
#     }


#     ## Write to NetCDF
#     compression = {"zlib": True, 
#                 "complevel": 1,
#                 "shuffle": True,}
#     encoding = {var: compression for var in ds.data_vars}
#     output_dir = Path(f"/discover/nobackup/projects/QEFM/data/rollout_outputs/{fmodel}/geos-fp-interp-no-mask/Y{yyyy}/M{mm}/D{dd}")
#     #output_dir = Path(f"/discover/nobackup/projects/QEFM/data/rollout_outputs/{fmodel}/0p25/Y{yyyy}/M{mm}/D{dd}")
#     output_dir.mkdir(parents=True, exist_ok=True)
#     fname = f"{fmodel}-prediction-geos_date-{tstamp}_res-1.0_levels-13_ens-mean.nc"
#     output_file = output_dir / fname
#     ds.to_netcdf(output_file, encoding=encoding, engine="netcdf4")



#     print("Finished \n", ds)
