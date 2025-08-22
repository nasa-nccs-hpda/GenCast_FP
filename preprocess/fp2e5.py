import os
import numpy as np
import xarray as xr
from datetime import datetime
import xesmf  as xe
from fp_to_era5 import *

if __name__ == "__main__":
    for day in range(1, 32):
        for hour in [0, 12]:
            dt = datetime(2024, 5, day, hour)
            Files =  discover_files(dt, outdir='./output/')
            print(Files)

            fp_Nx = xr.open_dataset(Files['fp_Nx'],engine='netcdf4')
            fp_Nv = xr.open_dataset(Files['fp_Nv'],engine='netcdf4')

            ai_Nx, ai_Np = fp_to_era5_xlevs ( fp_Nx, fp_Nv )

            ai_Ex = era5_dataset('GEOS-FP 2D Variables on ERA-5 Grid for AI/ML Modeling')
            regridder = xe.Regridder(ai_Nx, ai_Ex, "conservative") # do this once

            ai_Ex, ai_Ep = fp_to_era5_hgrid ( ai_Nx, ai_Np, regridder=regridder)

            ai_Ex.to_netcdf(Files['ai_Ex'],engine='netcdf4')
            ai_Ep.to_netcdf(Files['ai_Ep'],engine='netcdf4')    
