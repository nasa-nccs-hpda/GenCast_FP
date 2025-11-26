"""
Convert these collections from GEOS-FP

- asm_Nv
- asm_Nx

to pressure coordinates, renaming variables as ERA-5,
and converting to ERA-5 lat-lon grid.

Unlike standard GEOS-FP files, the output files are
extrapolated under the mountains using CDO.

TODO:
    - PEP8 compliance
"""

import os
import numpy  as np
import xarray as xr
import xesmf  as xe

from datetime import datetime

from gencast_fp.preprocess.eta2xprs import xEta2xprs
from gencast_fp.preprocess.eta2xprs import LOG

import gencast_fp.preprocess.eta as eta


# Pressure levels for AI/ML ERA-5 like file
# -----------------------------------------
e5_plevs  = np.array(
    [
        50, 100, 150, 200, 250, 300, 400,
        500, 600, 700, 850, 925, 1000
    ]
)  # hPa

GRAV = 9.80665

# Attribute fixer
# ---------------
e5_attrs = dict(
    w={'long_name': 'Vertical velocity', 'units': 'Pa s**-1',
        'standard_name': 'lagrangian_tendency_of_air_pressure'},
    q={'long_name': 'Specific humidity', 'units': 'kg kg**-1',
        'standard_name': 'specific_humidity'},
    t={'long_name': 'Temperature', 'units': 'K',
        'standard_name': 'air_temperature'},
    u={'long_name': 'U component of wind', 'units': 'm s**-1',
        'standard_name': 'eastward_wind'},
    v={'long_name': 'V component of wind', 'units': 'm s**-1',
        'standard_name': 'northward_wind'},
    p={'long_name': 'Mid-layer Pressure', 'units': 'Pa',
        'standard_name': 'air_pressure'},
    hgt={'long_name': 'Geopotential Height', 'units': 'm',
         'standard_name': 'geopotential_height'},
    z={'long_name': 'Geopotential', 'units': 'm**2 s**-2',
        'standard_name': 'geopotential'},
    sp={'long_name': 'Surface pressure', 'units': 'Pa',
        'standard_name': 'surface_air_pressure'},
    msl={'long_name': 'Mean sea level pressure', 'units': 'Pa',
         'standard_name': 'air_pressure_at_mean_sea_level'},
    t2m={'long_name': '2 metre temperature', 'units': 'K',
         'standard_name': 'air_temperature'},
    skt={'long_name': 'Skin temperature', 'units': 'K',
         'standard_name': 'surface_skin_temperature'},
    u10={'long_name': '10 metre U wind component', 'units': 'm s**-1',
         'standard_name': 'eastward_wind'},
    v10={'long_name': '10 metre V wind component', 'units': 'm s**-1',
         'standard_name': 'northward_wind'},
)


def _fixAttrs(v, ds, ref_attrs=e5_attrs):
    """
    Fix key attributes.
    """
    for a in ['long_name', 'units', 'standard_name']:
        if len(ds[v].shape) > 2:
            ds[v].attrs[a] = ref_attrs[v][a]


def _fixPoles_scalar(a5, ap):
    """
    Fix scalar variable at poles. Assumes first and last
    latitudinal points are +/- 90
    in both source (ap) and target (a5). Notice that ERA-5 has latitudinal
    grid north-to-south while FP is south-to-north.

    Indexing is assumed to be one of the following:
        - a(lev,lat,lon)
        - a(lat,lon)

    """
    shape = ap.shape
    if len(shape) == 3:    # 3D
        for k in range(shape[0]):
            a5[k, 0, :] = ap[k,-1, :].mean()  # north pole
            a5[k, -1, :] = ap[k, 0, :].mean()  # south pole
    elif len(shape) == 2:  # 2D
        a5[0, :]  = ap[-1, :].mean(axis=0)  # north pole
        a5[-1, :] = ap[0, :].mean(axis=0)  # south pole       
    else:
        print(ap.shape)
        raise ValueError('Invalid shape of input variable')


def _fixPoles_vector(u5, v5, up, vp):
    """
    Fix vector variables at poles. Assumes first and
    last latitudinal points are +/- 90
    in both source (ap) and target (a5). Notice that ERA-5 has latitudinal grid
    north-to-south while FP is south-to-north.

    This is implemented with a simple linear interpolation in longitude.

    """

    lon_p = up.lon
    lon_5 = u5.longitude

    shape = up.shape

    if len(shape) == 3:    # 3D
        for k in range(shape[0]):
            u5[k, 0].data[:] = np.interp(
                lon_5, lon_p, up[k,-1].data[:], period=360.)
            v5[k, 0].data[:] = np.interp(
                lon_5, lon_p, vp[k,-1].data[:], period=360.)        
            u5[k, -1].data[:] = np.interp(
                lon_5, lon_p, up[k, 0].data[:], period=360.)
            v5[k, -1].data[:] = np.interp(
                lon_5, lon_p, vp[k, 0].data[:], period=360.)

    elif len(shape) == 2:  # 2D
        u5[0].data[:] = np.interp(lon_5, lon_p, up[-1].data[:], period=360.)
        v5[0].data[:] = np.interp(lon_5, lon_p, vp[-1].data[:], period=360.)
        u5[-1].data[:] = np.interp(lon_5, lon_p, up[0].data[:], period=360.)
        v5[-1].data[:] = np.interp(lon_5, lon_p, vp[0].data[:], period=360.)

    else:
        print(up.shape)
        raise ValueError('Invalid shape of input variable')


def _scalar_vectors(ds):
    """
    Given a dataset, return list of scalar and vector variables.
    """

    # Identify U, V and scalars
    # -------------------------
    U, V, S = {}, {}, []
    for v in ds.data_vars:
        if len(ds[v].shape) < 3:
            continue
        std_name = ds[v].attrs['standard_name']
        if 'eastward' in std_name:
            U[std_name] = v
        elif 'northward' in std_name:
            V[std_name] = v
        else:
            S += [v,]

    # Vector pairs
    # ------------
    VP = []
    for u_ in U:
        v_ = u_.replace('eastward', 'northward')
        p = (U[u_], V[v_])
        VP += [p,]

    return (S, VP)


def _gat2s(template, time, expid='f5295'):
    """
    Expand GrADS style templates/
    """
    
    y4 = '%04d'%time.year
    m2 = '%02d'%time.month
    d2 = '%02d'%time.day
    h2 = '%02d'%time.hour
    n2 = '%02d'%time.minute
    
    return template.replace('%y4',y4).replace('%m2',m2).\
                    replace('%d2',d2).replace('%h2',h2).\
                    replace('%n2',n2).replace('%expid',expid)


def discover_files(time, outdir='./', expid='f5295'):
    """
    Return dictionary with FP file names on discover, given time.
        
    Input GEOS-FP Files:
        fp_Nv, fp_Nx, fp_Np : GEOS-FP horizontal grid, v=eta, p=pressure, x=surface
        
    Output AI/ML files: 
        ai_Np, ai_Nx    : GEOS-FP horizontal grid, p=pressure, x=surface
        ai_Ep, ai_Ex    : ERA-5 horizontal grid, p=pressure, x=surface

    Optional ERA-5 Files:
        e5_Ep, e5_Ex    : ERA-5 horizontal grid, p=pressure, x=surface
        
    """

    # ERA-5 and GEOS-FP directories on discover
    # -----------------------------------------
    t, x = time, expid
    sst_dir = '/gpfsm/dnb07/projects/p10/gmao_ops/fvInput/g5gcm/bcs/realtime/OSTIA_REYNOLDS/2880x1440/'
    e5_dirn = '/css/era5/static/'
    e5_dir2 = '/css/era5/surface_hourly/inst/Y%y4/M%m2/'
    fp_dirn = '/gpfsm/dnb06/projects/p174/%expid_fp/diag/Y%y4/M%m2/'

    F1 = dict(e5_Es = (e5_dirn+'era5_static-allvar.nc'),
              e5_Ex = _gat2s(e5_dir2+'era5_surface-inst_allvar_%y4%m2%d2_%h2z.nc',t,x),
              sst = _gat2s(sst_dir+'dataoceanfile_OSTIA_REYNOLDS_SST.2880x1440.%y4.data',t,x),
              fp_Nv = _gat2s(fp_dirn+'%expid_fp.inst3_3d_asm_Nv.%y4%m2%d2_%h2%n2z.nc4',t,x),
              fp_Nx = _gat2s(fp_dirn+'%expid_fp.inst3_2d_asm_Nx.%y4%m2%d2_%h2%n2z.nc4',t,x),
              )

    F2 = dict(ai_Nv = outdir+os.path.basename(F1['fp_Nv']).replace('asm_Nv','ai_Nv'),
              ai_Np = outdir+os.path.basename(F1['fp_Nv']).replace('asm_Nv','ai_Np'),
              ai_Nx = outdir+os.path.basename(F1['fp_Nx']).replace('asm_Nx','ai_Nx'),
              ai_Ex = outdir+os.path.basename(F1['fp_Nx']).replace('asm_Nx','ai_Ex'),
              ai_Ep = outdir+os.path.basename(F1['fp_Nv']).replace('asm_Nv','ai_Ev'),
             )

    F1.update(F2)
    return F1

def era5_dataset(title,e5_plevs=e5_plevs):
    """
    Return an empty dataset with ERA-5 coordinates and simple attributes.
    """

    # Horizontal coords hardwired
    # ---------------------------
    e5_lons = np.linspace(0, 359.75, num=1440)
    e5_lats = np.linspace(  90., -90., num=721)

    lon = xr.DataArray(e5_lons, dims='longitude',
                                attrs = dict ( long_name='Longitude',
                                               units='degrees_east',
                                               standard_name='longitude') )
    lat = xr.DataArray(e5_lats, dims='latitude',
                                attrs = dict ( long_name='Latitude',
                                               units='degrees_north',
                                               standard_name='latitude') )
    lev = xr.DataArray(e5_plevs, dims='pressure_level',
                                 attrs = dict ( long_name='Pressure',
                                                units='hPa',
                                                standard_name='air_pressure') )

    # Hardwire attributes
    # -------------------
    attrs = dict ( Title = title,
                   Conventions = 'CF-1.7',
                   Institution = 'Global Modeling and Assimilation Office, NASA/GSFC',
                   Contact = 'arlindo.m.dasilva@nasa.gov',
                   History = 'Created with fp_to_era5 on %s'%str(datetime.now())
                   )
                   
    
    coords = dict ( longitude=lon, latitude=lat, pressure_level=lev )

    return xr.Dataset(coords=coords, attrs=attrs)

def sst_dataset(title):
    """
    Return an empty dataset with SST coordinates and simple attributes.
    """

    # Horizontal coords hardwired
    # ---------------------------
    res = 0.125
    sst_lons = -179.935 + res * np.arange(2880)
    sst_lats = -89.9375 + res * np.arange(1440)

    lon = xr.DataArray(sst_lons, dims='longitude',
                                attrs = dict ( long_name='Longitude',
                                               units='degrees_east',
                                               standard_name='longitude') )
    lat = xr.DataArray(sst_lats, dims='latitude',
                                attrs = dict ( long_name='Latitude',
                                               units='degrees_north',
                                               standard_name='latitude') )

    # Hardwire attributes
    # -------------------
    attrs = dict ( Title = title,
                   Conventions = 'CF-1.7',
                   Institution = 'Global Modeling and Assimilation Office, NASA/GSFC',
                   History = 'Created with fp_to_era5 on %s'%str(datetime.now()),
    )

    coords = dict ( longitude=lon, latitude=lat)

    return xr.Dataset(coords=coords, attrs=attrs)

#........................................................................................
def fp_to_era5_xlevs ( fp_Nx, fp_Nv, plevs=e5_plevs ):
    
    """
    Given asm_Nv and asm_Nx collections from GEOS-FP given as xarray
    do this:

    - Restructure data sets to match ECMWF names, dropping variables
      not needed by GenCast.
    - Interpolate to pressure level, extrapolating to below the surface
      using the so-called "ECMWF method"

      ai_Nv, ai_Nx = fp_to_era5 ( fp_Nv, fp_Nx )

      On output, these datasets are still on GEOS-FP horizontal grid.
      
    """

    # Trim and rename 3D FP variables
    # -------------------------------
    #print('Original 3D Data Variables:\n ',list(fp_Nv.data_vars))
    fp_Nv = fp_Nv.drop_vars(['CLOUD', 'DELP', 'EPV', 'O3', 'QI', 'QL', 'QR', 'QS','RH']).\
                            rename(dict(OMEGA='w', PHIS='zs', SLP='msl', PS='sp', 
                                        H='hgt', PL='p', QV='q', T='t', U='u', V='v'))

    # Trim and rename GEOS-FP 2D Variables
    # ------------------------------------
    #print('Original 2D Data Variables:\n ',list(fp_Nx.data_vars))
    fp_Nx = fp_Nx.drop_vars(['DISPH', 'HOURNORAIN', 'PS', 'QV10M', 'QV2M', 'SLP', 'T10M',
                         'T2MMAX', 'T2MMIN', 'TO3', 'TOX','TQI', 'TQL', 
                         'TQV', 'TROPPB', 'TROPPT', 'TROPPV', 'TROPQ', 'TROPT',
                         'U2M', 'U50M', 'V2M', 'V50M']).\
                         rename(dict( TPRECMAX='tp', T2M='t2m', TS='skt',
                                      U10M='u10', V10M='v10'))
    
    #print('Kept 2D Data Variables:\n ',list(fp_Nx.data_vars))
    
    # Fix 3D variable attributes
    # --------------------------
    fp_Nv['hgt'].attrs['standard_name'] = 'geopotential_height'
    fp_Nv['p'].attrs['standard_name'] = 'air_pressure'
    for v in [ 'w', 'q', 't', 'u', 'v' ]:
        _fixAttrs ( v, fp_Nv)
    
    # Fix 2D variable atributes in 3D file
    # --------------------------
    fp_Nv['zs'].attrs['long_name'] = 'Surface Geopotential'
    fp_Nv['zs'].attrs['standard_name'] = 'surface_geopotential'
    fp_Nv['zs'].attrs['units'] = 'm2 s-2'
    for v in [ 'sp', 'msl' ]:
        _fixAttrs ( v, fp_Nv)
    
    # Fix 2D Variable attribute in 2D file
    # ------------------------------------
    fp_Nx['tp'].attrs['standard_name'] = 'precipitation_amount'
    for v in  ['t2m', 'skt', 'u10', 'v10']:
        _fixAttrs ( v, fp_Nx)

    # Move 2D vars from 3D file
    # -------------------------
    drop_list = []
    for v in fp_Nv.data_vars:
        if len(fp_Nv[v].shape) == 3:
            fp_Nx[v] = fp_Nv[v]
            drop_list += [v,]
    drop_list.remove('sp')  # needed by CF to specify vertical coordinate
    
    #drop_list.remove('zs')  # interpolation needs this in the 3D file
    #drop_list.remove('hyam')  # interpolation needs this in the 3D file
    #drop_list.remove('hybm')  # interpolation needs this in the 3D file
    
    fp_Nv = fp_Nv.drop_vars(drop_list)

    # CF compliant eta coordinate encoding
    # ------------------------------------
    ak_, bk_ = np.array(eta.ak['72']), np.array(eta.bk['72']) # edge levels
    ak = ( ak_[:-1] + ak_[1:] ) / 2.
    bk = ( bk_[:-1] + bk_[1:] ) / 2.
    
    fp_Nv['lev'] = bk + ak/101325.
    fp_Nv.lev.attrs['units'] = '1'
    fp_Nv.lev.attrs['standard_name'] = 'atmosphere_hybrid_pressure_coordinate'
    fp_Nv.lev.attrs['formula_terms'] = 'ap: ak b: bk ps: sp'
    fp_Nv['ak'] = xr.DataArray(ak,dims=['lev'], attrs={'units':'pa'} )
    fp_Nv['bk'] = xr.DataArray(bk,dims=['lev'], attrs={'units':'1'} )
    
    # Prepare pressure coordinate dataset
    # -----------------------------------
    fp_Np = xEta2xprs(fp_Nv, fp_Nx, plevs, method=LOG)

    # Compute geopotential from heights
    # ---------------------------------
    fp_Np['z'] = GRAV * fp_Np['hgt']
    
    # Fix 3D Variable attribute in 2D file
    # ------------------------------------
    for v in fp_Np.data_vars:
        _fixAttrs ( v, fp_Np )

    # Return sanitized 2D and 3D datasets
    # -----------------------------------
    return fp_Nx, fp_Np

#--
def fp_to_era5_hgrid ( ai_Nx, ai_Np,
                       regridder=None, method="conservative"):
    """
    Regrid from GEOS-FP to ERA-5 horizontal grid using xESMF.

    ai_Ex, ai_Ep = fp_to_era5_hgrid ( ai_Nx, ai_Np, regridder)
    
    For efficiency, the regridder can be computed once and provided on
    input.

    Note: poles are specialed handled.
    
    """

    # Create regridder if not provided
    # --------------------------------
    if regridder is None:
        ai_Ex = era5_dataset('GEOS-FP 2D Variables on ERA-5 Grid for AI/ML Modeling')
        regridder = xe.Regridder(ai_Nx, ai_Ex, method)

    # Regrid
    # ------
    ai_Ex = regridder(ai_Nx,keep_attrs=True) 
    ai_Ep = regridder(ai_Np,keep_attrs=True)

    # Fix poles: 2D
    # -------------
    nt = ai_Ex.time.size
    S, V = _scalar_vectors(ai_Ex)
    for t in range(nt):
        for s in S:
            _fixPoles_scalar( ai_Ex[s][t], ai_Nx[s][t])
        for u, v in V:
            _fixPoles_vector( ai_Ex[u][t], ai_Ex[v][t],
                              ai_Nx[u][t], ai_Nx[v][t])
    # Fix poles: 3D
    # -------------
    nt = ai_Ex.time.size
    S, V = _scalar_vectors(ai_Ep)
    for t in range(nt):
        for s in S:
            _fixPoles_scalar( ai_Ep[s][t], ai_Np[s][t])
        for u, v in V:
            _fixPoles_vector( ai_Ep[u][t], ai_Ep[v][t],
                              ai_Np[u][t], ai_Np[v][t])

    return (ai_Ex, ai_Ep)
            
#....................................................................................
if __name__ == "__main__":


    # Directories
    # -----------
    e5_dir3d = '/css/era5/pressure_hourly/inst'
    e5_dir2d = '/css/era5/surface_hourly/inst'
    fp_dirn = '/gpfsm/dnb06/projects/p174/f5295_fp/diag'

    e5_fileNv = e5_dir3d + '/Y%y4/M%m2/era5_atmos-inst_allvar_20241212_12z.nc'
    e5_fileNx = e5_dir2d + '/Y%y4/M%m2/era5_surface-inst_allvar_20241212_12z.nc'

    fp_fileNv = fp_dirn +  '/Y2024//M12/f5295_fp.inst3_3d_asm_Nv.%y4%m2%d2_%h2%n2z.nc4'
    fp_fileNx = fp_dirn +  '/Y2024//M12/f5295_fp.inst3_2d_asm_Nx.%y4%m2%d2_%h2%n2z.nc4'

    aiml_Ex = os.path.basename(fp_fileNx).replace('asm_Nx','aiml_Ex')
    aiml_Ep = os.path.basename(fp_fileNv).replace('asm_Nv','aiml_Ev')

    # Lazy load FP files
    # ------------------
    fp_Nx = xr.open_dataset(fp_fileNx,engine='netcdf4')
    fp_Nv = xr.open_dataset(fp_fileNv,engine='netcdf4')


    
    # Save as NetCDF files
    # --------------------
    # if save_intermediate:
    #     aiml_Nv = os.path.basename(fp_fileNv).replace('asm_Nv','aiml_Nv')
    #     fp_Nv.to_netcdf(aiml_Nv,engine='netcdf4')
    #     aiml_Nx = os.path.basename(fp_fileNx).replace('asm_Nx','aiml_Nx')
    #     fp_Nx.to_netcdf(aiml_Nx,engine='netcdf4')


