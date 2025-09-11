"""
Simple wrapper around eta2prs Fortran extension.

Arlindo da Silva, April 2025.

"""

import numpy  as np
import xarray as xr
import eta
import eta2xprs_   # f2py extension

HEIGHT, TEMPERATURE, OTHER = -1, 1, 0
LINEAR, LOG, LOGLOG = 1, 2, 3

UNDEF = 1.0E10

def xEta2xprs(fp_Nv, fp_Nx, plevs, method, 
              ti=0, eta_coeffs=None, extrapolate=True, verbose=False):
    """
    Convert all data variables in the xarray dataset *fp_Nv* from eta to pressure coordinates.

    fp_Np = xEta2prs(fpNv, fp_Nx, plevs, method, ...)

    On Input
    --------
    fp_Nv          xr.Dataset      3d variables
    fp_Nx          xr.Dataset      2d variables

    
    plevs:         (plev)          Pressure levels to interpolate to [hPa]

    method         Method dor vertical interpolation
                      1  LINEAR
                      2  LOG
                      3  LOG LOG
                      
    vflag:         Variable flag: 
                     -1  means geopotential height (Z)
                     +1  means temperature (T)
                      0  any other variable

    eta_coeffs     Eta coefficients [ak, bk]. Dimension must match
                   the variable vertical dimension. If the 3D variables are
                   mid-layer, the ak/bk must be mid-layer. If the 3D variables are
                   edge-level variables, then ak/bk must be edge-level.
    
    extrapolate:   Extrapolation flag: 
                      False - do not extrapolate
                      True  - extrapolate using ECMWF method

    ti             int  time index
                      

    On Output:
    ----------
    f_prs          (plev,lat,lon)  array interpolated to pressure

    """

    # Coordinates for pressure coordinate dataset
    # -------------------------------------------
    prs = xr.DataArray(plevs, dims='pressure_level',
                              attrs = dict ( long_name='Pressure',
                                             units='hPa',
                                             standard_name='air_pressure') )
    p_coords = fp_Nv.coords.copy()
    p_coords['lev'] = prs  # will rename later

    # Create p-coord dataset
    # ----------------------
    Title, Filename = fp_Nv.attrs['Title'], fp_Nv.attrs['Filename']
    attrs = fp_Nv.attrs.copy()
    attrs['Title'] = Title.replace('Model-Level','Pressure-Level')
    attrs['Filename'] = Filename.replace('asm_Nv','asm_Np')

    fp_Np = xr.Dataset(coords=p_coords, attrs=attrs)

    # Reuse eta coefficients
    # ----------------------
    km = len(fp_Nv.coords['lev'])
    if eta_coeffs is None:
        try:
            ak, bk = np.array(eta.ak[str(km)]), np.array(eta.bk[str(km)])
            ak_, bk_ = (ak[0:-1]+ak[1:])/2., (bk[0:-1]+bk[1:])/2.             # mid-layer
        except:
            ak_, bk_ = np.array(eta.ak[str(km-1)]), np.array(eta.bk[str(km-1)]) # edges
        eta_coeffs = [ak_, bk_]
    
    # For every data variables
    # ------------------------
    for v in fp_Nv.data_vars:

        V = fp_Nv[v]
        
        # Special 2D variables
        # --------------------
        if len(V.shape) == 1:
            continue  # these are likely ak, bk
        elif len(V.shape) < 4:
            if verbose:
                print('[ ] Skipping eta to prs for 2D variable <%s>'%v)
            fp_Np[v] = V
            
        # 3D vars
        # -------
        else:

            if verbose:
                print('[x] Doing eta to prs for 3D variable <%s>'%v)

            # Use standard name to deduce vflag
            # ---------------------------------
            if V.attrs['standard_name'] == 'air_temperature':
                vflag = TEMPERATURE
            elif 'geopotential' in V.attrs['standard_name']:
                vflag = HEIGHT
            else:
                vflag = OTHER
            
            v_data = Eta2xprs(v, fp_Nv, fp_Nx, plevs, method, vflag,
                              ti, eta_coeffs, extrapolate)
            
            fp_Np[v] = xr.DataArray(v_data, dims=['time','lev','lat', 'lon'])
            
    fp_Np = fp_Np.rename({'lev':'pressure_level'}) # to be clear
    
    return fp_Np


def Eta2xprs(v, fp_Nv, fp_Nx, plevs, method, vflag,
             ti=0, eta_coeffs=None, extrapolate=True):
    """
    Interpolate from eta (hybrid sigma pressure) to pressure levels
    using the so-called ECMWF method. Xarray interface.

    f_prs = Eta2prs(v, fp_Nv, fp_Nv, plevs, method, vflag, extrapolate=True)

    All arrays assumed to be numpy arrays.
    
    On Input
    --------
    v              str             variable name
    fp_Nv          xr.Dataset      3d variables
    fp_Nx          xr.Dataset      2d variables

    
    plevs:         (plev)          Pressure levels to interpolate to [hPa]

    method         Method dor vertical interpolation
                      1 - LINEAR
                      2 - LOG
                      3 - LOG LOG
                      
    vflag:         Variable flag: 
                     -1  means geopotential height (Z)
                     +1  means temperature (T)
                      0  any other variable

    eta_coeffs     Eta coefficients [ak, bk]. Dimension must match
                   the variable *v*. If *v* is mid-layer, the ak/bk
                   must be mid-layer. If *v* is an edge-level variable, then
                   ak/bk must be edge-level.

    extrapolate:   Extrapolation flag: 
                      False - do not extrapolate
                      True  - extrapolate using ECMWF method

    ti             int  time index
                      

    On Output:
    ----------
    f_prs          (plev,lat,lon)  array interpolated to pressure

    
    """

    # Dimensions
    # ----------
    tm, im, jm, km = fp_Nv[v].shape
    kp = len(plevs)

    f, ps, ts, phis = (fp_Nv[v][ti].values.T, fp_Nx['sp'][ti].values.T, 
                       fp_Nx['skt'][ti].values.T, fp_Nx['zs'][ti].values.T)

    # Mid-layer or edge levels
    # ------------------------
    if eta_coeffs is None:
        try:
            ak, bk = np.array(eta.ak[str(km)]), np.array(eta.bk[str(km)])
            ak_, bk_ = (ak[0:-1]+ak[1:])/2., (bk[0:-1]+bk[1:])/2.             # mid-layer
        except:
            ak_, bk_ = np.array(eta.ak[str(km-1)]), np.array(eta.bk[str(km-1)]) # edges
    else:
        ak_, bk_ = eta_coeffs
        
    undef = 1.0E10 # undefined value

    # a_prs = eta2prs(a_eta, ak, bk, ps, ts, phis, plevs,
    #                 xflag, vflag, method, undef)
    
    f_p = eta2xprs_.eta2xprs ( f, ak_, bk_, ps, ts, phis, plevs, 
                              extrapolate, vflag, method, undef)
    V = f_p.T
    V = V.reshape((1,)+V.shape) # restore time dimensions
    
    return V

def eta2xprs(f_eta,ps,ts,phis,plevs,method,vflag,
             eta_coeffs=None, extrapolate=True):
    """
    Interpolate from eta (hybrid sigma pressure) to pressure levels
    using the so-called ECMWF method. Numpy array interface.

    f_prs = eta2prs(f_eta,ps,ts,phis,plevs,method,vflag,extrapolate=True)

    All arrays assumed to be numpy arrays.
    
    On Input
    --------
    f_eta          (lev,lat,lon)   numpy array to be interpolated to pressure
    ps             (lat,lon)       surface pressure [Pa]
    ts             (lat,lon)       near ground temperature [K]
    phis           (lat,lon)       geopotential height [m2 s-2]
    
    plevs:         (plev)          Pressure levels to interpolate to [hPa]

    method         Method dor vertical interpolation
                      1 - LINEAR
                      2 - LOG
                      3 - LOG LOG
                      
    extrapolate:   Extrapolation flag: 
                      False - do not extrapolate
                      True  - extrapolate using ECMWF method
                      
    vflag:         Variable flag: 
                     -1  means geopotential height (Z)
                     +1  means temperature (T)
                      0  any other variable

    eta_coeffs     Eta coefficients [ak, bk]. Dimension must match
                   verticla dimension of *f_eta*. If *f_eta* is mid-layer, 
                   the ak/bk must be mid-layer. If *f_eta* is an edge-level variable, 
                   then ak/bk must be edge-level.

    On Output:
    ----------
    f_prs          (plev,lat,lon)  numpy array interpolated to pressure

    
    """

    # Mid-layer or edge levels
    # ------------------------
    if eta_coeffs is None:
        try:
            ak, bk = np.array(eta.ak[str(km)]), np.array(eta.bk[str(km)])
            ak_, bk_ = (ak[0:-1]+ak[1:])/2., (bk[0:-1]+bk[1:])/2.               # mid-layer
        except:
            ak_, bk_ = np.array(eta.ak[str(km-1)]), np.array(eta.bk[str(km-1)]) # edges
    else:
        ak_, bk_ = eta_coeffs
        
    undef = 1.0E10 # undefined value

    # a_prs = eta2prs(a_eta, ak, bk, ps, ts, phis, plevs,
    #                 xflag, vflag, method, undef)
    
    V = eta2xprs_.eta2xprs ( f.T, ak_, bk_, ps.T, ts.T, phis.T, plevs, 
                             extrapolate, vflag, method, UNDEF)

    return V.T

