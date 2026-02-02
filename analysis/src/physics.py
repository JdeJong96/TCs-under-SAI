#!/usr/bin/env python

import numpy as np
import xarray as xr

A = 6371000 # radius earth [m]
PI = np.pi  # pi
G = 9.81    # acceleration due to gravity [m/s^2]


def cumdiff(da, dim):
    """calculate difference along dimension keeping the first element
    a cumulative sum on the result will reproduce the original DataArray
    

    da : xr.DataArray
        input DataArray
    dim : str
        dimension

    returns:
        xr.DataArray : resulting difference array (same shape as da)
    """

    return xr.concat((
        da.isel({dim: 0}),
        da.diff(dim, label='upper')
    ), dim=dim)


def meridional_mass_streamfunction(v, dp, lat):
    """meridional mass streamfunction calculation (zonal mean)
    If longitude is a dimension of v, the result should be averaged afterwards.

    v : xr.DataArray
        meridional velocity [m/s]
    dp : xr.DataArray
        vertical pressure difference [Pa]
    lat : xr.DataArray
        latitude [deg N]

    returns:
    xr.DataArray : (same shape as v)
        meridional mass streamfunction [kg/s]  
    """
    result = 2*PI*A/G * np.cos(np.deg2rad(lat)) * (v * dp).cumsum('lev')
    result.attrs.update({
        'long_name':'meridional mass streamfunction',
        'units': 'kg/s'
    })
    result.name = 'MMF'
    return result.astype(v.dtype)


def windshear_250_850(ds):
    """Calculate vertical windshear between 250 and 850 hPa

    input:
    ds : xr.Dataset
        dataset containing zonal and meridional wind at 250 and 850 hPa

    returns:
    xr.DataArray:
        absolute vertical wind difference between level 1 and 2
    """
    
    vws = np.sqrt((ds.U250 - ds.U850)**2 + (ds.V250 - ds.V850)**2)
    vws.name = 'VWS'
    vws = vws.assign_attrs({
        'units': ds.U250.units,
        'long_name': 'vertical wind shear between 250 and 850 hPa'
    })

    return vws