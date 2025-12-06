#!/usr/bin/env python

import numpy as np
import xarray as xr
from .load_SAIdata import Cases

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


def meridional_mass_flux(V, P, lat, a=6371000, g=9.81):
    """meridional mass flux calculation

    V : xr.DataArray
        meridional velocity [m/s]
    P : xr.DataArray
        pressure [Pa]
    lat : xr.DataArray
        latitude [deg N]
    a : int | float
        radius of earth [m]
    g : int | float
        acceleration due to gravity [m/s^2]

    returns:
        xr.DataArray : meridional mass flux (same shape as V)
    """

    dP = cumdiff(P, 'lev')
    MMF = 2*np.pi*a/g * np.cos(np.deg2rad(lat)) * (V * dP).cumsum('lev')

    return MMF


if __name__ == '__main__':
    ds = Cases('hres.ref.1').select('atm', 'h0').open_mfdataset(decode_times=False, chunks={})

    #ds = ds.isel(time=0)

    P = (ds.hyam * ds.P0 + ds.hybm * ds.PS)

    MMF = meridional_mass_flux(ds.V, P, ds.lat)


