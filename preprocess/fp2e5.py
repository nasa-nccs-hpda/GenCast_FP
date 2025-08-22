## this is the driver code for converting GEOS-FP data to ERA5 format
import os
import numpy as np
import xarray as xr
from datetime import datetime
import xesmf  as xe
from fp_to_era5 import *
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert GEOS-FP data to ERA5 format.')
    parser.add_argument('--year', type=int, required=True, help='Year of the data to process')
    parser.add_argument('--month', type=int, required=True, help='Month of the data to process')
    parser.add_argument('--day', type=int, required=True, help='Day of the data to process')
    
    parser.add_argument('--outdir', type=str, default='./output/', help='Output directory for the converted files')
    parser.add_argument('--expid', type=str, default='fp2e5', help='Experiment ID for the output files')

    args = parser.parse_args()
    yyyy = args.year
    mm  = args.month
    dd = args.day
    outdir = args.outdir
    expid = args.expid

    regridder = None
    for hour in [0, 12]:
        dt = datetime(yyyy, mm, dd, hour)
        Files = discover_files(dt, outdir=outdir, expid=expid)
        print(Files)

        fp_Nx = xr.open_dataset(Files['fp_Nx'], engine='netcdf4')
        fp_Nv = xr.open_dataset(Files['fp_Nv'], engine='netcdf4')

        ai_Nx, ai_Np = fp_to_era5_xlevs(fp_Nx, fp_Nv)

        ai_Ex = era5_dataset('GEOS-FP 2D Variables on ERA-5 Grid for AI/ML Modeling')
        if not regridder:
            regridder = xe.Regridder(ai_Nx, ai_Ex, "conservative")
        
        ai_Ex, ai_Ep = fp_to_era5_hgrid(ai_Nx, ai_Np, regridder=regridder)

        ai_Ex.to_netcdf(Files['ai_Ex'],engine='netcdf4')
        ai_Ep.to_netcdf(Files['ai_Ep'],engine='netcdf4')  



if __name__ == "__main__":
    main()
   
