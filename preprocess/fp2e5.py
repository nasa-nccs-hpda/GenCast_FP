## this is the driver code for converting GEOS-FP data to ERA5 format
import os
import numpy as np
import xarray as xr
import pandas as pd
import xesmf  as xe
from fp_to_era5 import *
import argparse

def main():
    parser = argparse.ArgumentParser(description='Convert GEOS-FP data to ERA5 format.')
    parser.add_argument('--start_date', type=str, required=True, help='Start date to process (YYYY-MM-DD)')
    parser.add_argument('--end_date', type=str, required=True, help='End date to process (YYYY-MM-DD)')
    
    parser.add_argument('--outdir', type=str, default='./output/', help='Output directory for the converted files')
    parser.add_argument('--expid', type=str, default='f5295', help='Experiment ID for the output files')

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    dates = pd.date_range(start=args.start_date, end=args.end_date, freq='12H')


    outdir = args.outdir
    expid = args.expid

    regridder = None
    for _, dates in dates.groupby(dates.date).items():
        daily_Ex = []
        daily_Ep = []
        for dt in dates:
            Files = discover_files(dt, outdir=outdir, expid=expid)
            print(Files)

            fp_Nx = xr.open_dataset(Files['fp_Nx'], engine='netcdf4')
            fp_Nv = xr.open_dataset(Files['fp_Nv'], engine='netcdf4')

            ai_Nx, ai_Np = fp_to_era5_xlevs(fp_Nx, fp_Nv)

            ai_Ex = era5_dataset('GEOS-FP 2D Variables on ERA-5 Grid for AI/ML Modeling')
            if not regridder:
                regridder = xe.Regridder(ai_Nx, ai_Ex, "conservative")
            
            ai_Ex, ai_Ep = fp_to_era5_hgrid(ai_Nx, ai_Np, regridder=regridder)

            # add time dimension
            ai_Ex = ai_Ex.expand_dims(time=[dt])
            ai_Ep = ai_Ep.expand_dims(time=[dt])

            daily_Ex.append(ai_Ex)
            daily_Ep.append(ai_Ep)
        

        # concat along time
        ai_Ex_day = xr.concat(daily_Ex, dim="time")
        ai_Ep_day = xr.concat(daily_Ep, dim="time")

        # merge into single dataset for the day
        ai_day = xr.merge([ai_Ex_day, ai_Ep_day])
        print((ai_day))


if __name__ == "__main__":
    main()
   
