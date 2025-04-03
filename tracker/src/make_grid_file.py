#!/usr/bin/env python

import os
import numpy as np
import xarray as xr

file_05deg = "/projects/0/nwo2021025/archive/mres_b.e10.B2000_CAM5.f05_t12.001/atm/hist/mres_b.e10.B2000_CAM5.f05_t12.001.cam2.h0.2050-01.nc"
savepath = "/home/jasperdj/files_rene/Atmosphere_0_5_DX_DY_AREA.nc"
Rearth = 6371000 # Radius earth [m]

ds = xr.open_dataset(file_05deg)

lat = ds.lat.values
lon = ds.lon.values
gw = ds.gw.values

dlat = [(lat[0]+lat[1])/2-(-90), *[(lat[i+1]-lat[i-1])/2 for i in range(1,len(lat)-1)], 90-(lat[-1]+lat[-2])/2]
dlon = [(lon[(i+1)%len(lon)]-lon[i-1])/2%180 for i in range(len(lon))]

dx = np.cos(np.deg2rad(lat))[:,None] * np.deg2rad(dlon)[None,:] * Rearth
dy = np.ones_like(dx) * (np.deg2rad(dlat)*Rearth)[:,None]
area = np.ones_like(dx) * gw[:,None]/gw.sum()*4*np.pi*Rearth**2*dlon/360

ds.close()

ds_grid = xr.Dataset(
    data_vars={
        'DX': (['lat','lon'], dx, {'long_name':'Zonal length of grid cell','units':'m'}),
        'DY': (['lat','lon'], dy, {'long_name':'Meridional length of grid cell','units':'m'}),
        'AREA': (['lat','lon'], area, {'long_name':'Area of grid cell', 'units':'m^2'}),
    }, 
    coords={
        'lat': ('lat', lat, {'long_name':'Array of latitudes', 'units':'Degrees N'}),
        'lon': ('lon', lon, {'long_name':'Array of longitudes', 'units':'Degrees E'}),
    }, 
    attrs={
        'description':'Properties of the CESM 0.5 degree atmospheric grid'
    }
)

if os.path.exists(savepath):
    print(f"{savepath} already exists, delete first to create new file")
else:
    print(f"writing to {savepath}")
    ds_grid.to_netcdf(savepath)
