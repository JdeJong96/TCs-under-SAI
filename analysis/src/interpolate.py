#!/usr/bin/python
"""Vertical interpolation from atmospheric hybrid to pressure levels.

This script will read data on hybrid levels and interpolate them linearly onto 
pressure levels in log-pressure coordinates. Pressure is treated separately.
"""
import os
# if CPU usage does not exceed 100% (e.g. ~1500% for 16 cores), try:
#os.environ["OMP_NUM_THREADS"] = os.environ["SLURM_CPUS_PER_TASK"]
import sys
import glob
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from numba import float32, float64, guvectorize

# pressure levels [Pa] with similar resolution to original data
# see also analysis/notebooks/pres_levels.ipynb
PLEV = np.array([
    364.346569404006, 759.481964632869, 1435.6632251292501, 2461.22200042009, 
    3826.82997733355, 5459.54797416925, 7201.24505460262, 8782.123029232029, 
    10331.712663173701, 12154.7240763903, 14299.403876066199, 16822.5079774857, 
    20000.0, 22500.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 50000.0, 
    55000.0, 60000.0, 65000.0, 70000.0, 75000.0, 80000.0, 85000.0, 87500.0, 
    90000.0, 92500.0, 95000.0, 97000.0, 98500.0, 100000.0
])


def pressure_from_hybrid(ds, cell_interface=False):
    """Calculate atmospheric pressure from hybrid coefficients

    ds : xr.Dataset
        dataset from which to obtain hybrid coefficients, 
        reference and surface pressure
    cell_interface : Bool
        if True, calculate pressure at vertical cell interfaces, 
        else use midpoints (default)

    returns:
    xr.DataArray : 3D atmospheric pressure, same units as ds.PS
    """

    if cell_interface:
        P = ds.hyai * ds.P0 + ds.hybi * ds.PS
    else:
        P = ds.hyam * ds.P0 + ds.hybm * ds.PS
    P = P.assign_attrs(ds.PS.attrs)
    P.attrs.update({'long_name':'Pressure'})
    
    return P


@guvectorize(
    "(float64[:], float64[:], float64[:], float32[:])",
    " (n), (n), (m) -> (m)",
    nopython=True, target='parallel'
)
def interp1d_gu(f, x, xi, out):
    """Interpolate one-dimensional field f(x) to xi in ln(x) coordinates."""
    i, imax, x0, f0 = 0, len(xi), x[0], f[0]
    while (xi[i] < x0) and (i < imax):
        out[i] = np.nan      
        i = i + 1 
    for x1,f1 in zip(x[1:], f[1:]):
        while (xi[i] <= x1) and (i < imax):
            out[i] = (f1-f0)/np.log(x1/x0)*np.log(xi[i]/x0)+f0
            i = i + 1
        x0, f0 = x1, f1
    while i < imax:
        out[i] = np.nan
        i = i + 1


def interpolate(da, p, pi, dim, chunks={}):
    """Apply 1D interpolation function on xarray objects

    input:
    da : xr.DataArray
        field to interpolate
    p : xr.DataArray
        pressure
    pi : ArrayLike
        new pressure levels (1D)
    dim : Str
        dimension used for interpolation
    chunks : dict
        passed to da.chunk. if empty, existing chunks are retained

    returns:
    xr.DataArray
        da interpolated to new pressure levels
        dim is replaced by 'plev'
    """
    pi = xr.DataArray(name='plev', data=sorted(pi), dims='plev')
    dims = list(da.dims)
    dims[dims.index(dim)] = pi.dims[0]
    da = da.chunk(chunks)
    p = p.chunk(chunks)
    
    interped = xr.apply_ufunc(
        interp1d_gu, da, p, pi,
        input_core_dims=[[dim],[dim],pi.dims],
        output_core_dims=[pi.dims],
        exclude_dims=set((dim,)),
        dask='parallelized',
        on_missing_core_dim='drop'
    ).transpose(*dims).assign_coords({pi.name:pi}).assign_attrs(da.attrs)
    interped.name = da.name

    return interped


@guvectorize(
    "(float32[:], float64[:], float64, float64, float64[:])",
    " (n), (m), (), () -> ()",
    nopython=True, target='parallel'
)
def integrate1d_gu(f, pi, p_min, p_max, out):
    """Vertically integrate field f from p_min to p_max in pressure coordinates."""
    out[0] = 0.0
    nlev = len(f)
    for k in range(nlev):
        p_bnd = pi[k], pi[k+1]
        dp = min(p_bnd[1],p_max) - max(p_bnd[0],p_min)
        if dp <= 0:
            continue
        out[0] += f[k] * dp


def integrate(da, pi, p_min, p_max, da_dim='lev', pi_dim='ilev'):
    """Apply 1D integration function on xarray objects

    input:
    da : xr.DataArray[float32], shape lev
        field to integrate
    pi : xr.DataArray[float64], shape lev + 1
        pressure at vertical cell interfaces
        (e.g. pi = pressure_from_hybrid(ds, cell_interface=True))
    p_min : float64
        lower integration limit of pressure
    p_max : float64
        upper integration limit of pressure
    da_dim : Str
        vertical dimension in da
    pi_dim : Str
        vertical dimension in pi

    returns:
    xr.DataArray[float32]
        integral of da dp between p_min and p_max.
    """
    return xr.apply_ufunc(
        integrate1d_gu, da, pi, p_min, p_max,
        input_core_dims=[[da_dim],[pi_dim],[],[]],
        output_core_dims=[[]],
        dask='parallelized',
        on_missing_core_dim='drop'
    )


def main(ds):
    """Divide data into chunks, interpolate and average"""
    ds = ds.stack(x=('ens','time'))

if __name__ == '__main__':
	main()
