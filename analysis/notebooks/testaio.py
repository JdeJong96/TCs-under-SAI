import ast
import asyncio
import os
import glob
import numpy as np
import datetime
import xarray as xr
import nc_time_axis
import itertools
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
import shapely.geometry as sgeom
import cartopy.crs as ccrs
import cftime
from scipy.interpolate import interp1d
from src.load_SAIdata import Cases


async def add_data(arr, index, value):
    print(f"adding {value} to id {index}")
    await asyncio.sleep(1)
    arr[index] += value
    print(f"new array: {arr}")
    return


def preprocess_precipitation(data):
    """convert to 6-hourly average precipitation in mm/hour and center time stamps"""
    da = data.PRECT.copy()
    da.data *= 3.6e6
    da.attrs['units'] = 'mm/hour'
    dt = (da.time[1] - da.time[0]).dt.total_seconds().item() # time step (s)
    if 'time: mean' in getattr(da,'cell_methods',''): # set time to center of interval
        da['time'] = ('time', (da.time - datetime.timedelta(seconds=dt/2)).data, da.time.attrs)
    if dt == 10800: # convert 3hrly to 6hrly average
        t = da.time.coarsen(time=2).mean(keep_attrs=True) 
        da = da.coarsen(time=2).mean(keep_attrs=True)
        da['time'] = t
    assert ((da.time[1]-da.time[0]).dt.total_seconds()==21600, 
            f'expected 6-hourly CAM pcip... {da.time[:3].values=}')
    return da


def open_track_dataset(datadir, fix_gaps=False):
    """read in and merge data from all TC tracks

    datadir : str
        folder where all tracker data are present
        JdJ: ['../../tracker/data/','/home/jasperdj/files_rene/']
    """
    open_kwargs = {
        'Reference': ('RCP',2002,'ref'), 
        'RCP8.5': ('RCP',2092,'rcp'), 
        'SAI2050': ('SAI',2092,'sai')
    }
    
    ds = {exp: [] for exp in open_kwargs}
    for exp, (name, year, tag) in open_kwargs.items():
        maxtn = 0
        num_days = 0
        for n in range(1,6):
            fname = f'{name}.started_{year}.00{n}/TC_tracker_results.{tag}.started_{year}.00{n}.nc'
            fname = os.path.join(datadir, fname)
            if os.path.exists(fname):
                dsi = xr.open_dataset(fname, decode_cf=False)
                dsi = (dsi.assign_coords(id=('id',dsi.id.data+maxtn,dsi.id.attrs))
                       .assign_coords(ens=('id',np.ones(dsi.id.size)*n,{'long_name':'ensemble member'})))
                ds[exp].append(dsi)
                maxtn = dsi.id.max().item()
                times = dsi.TC_tracks.isel(data=0)
                num_days += dsi.num_days
                print(f'{tag.upper()}.00{n}: {dsi.num_days/365:.3f} years, {len(dsi.id)} tracks')
        ds[exp] = xr.concat(ds[exp], data_vars='minimal', dim='id')
        ds[exp]['num_days'] = num_days.assign_attrs({'long_name':'total number of analysed days'})
    
    
    for k,v in ds.items():
        for i,(desc) in enumerate(v.data.description):
            shortname = desc[:desc.index(':')]
            v[shortname] = v.TC_tracks.isel(data=i).assign_attrs(
                ast.literal_eval(desc[desc.index(':')+1:]))
            v[shortname].data[v[shortname]>1e30] = np.nan
        time = v.time.copy()
        tmask = np.isnan(time)
        time.data = cftime.num2date(time.fillna(0), time.units, time.calendar)
        ctime = time - datetime.timedelta(hours=1, minutes=30)
        v['year'] = ctime.dt.year.where(~tmask, np.nan) # at center of CESM time bounds
        v['cftime'] = time.where(~tmask, np.nan) #.where(~tmask, np.nan) + 
        v['PRECT'].data = v['PRECT']*3.6e6 # m/s to mm/hour
        v['PRECT'].attrs.update({'units':'mm/hour'})

    if not fix_gaps:
        return ds
    
    # There are a handful of NaNs in longitude in otherwise fine trajectories
    # fixing manually...
    for exp in ds:
        for tid in ds[exp].id.astype('int').data:
            dsi = ds[exp].sel(id=tid)
            if np.isnan(dsi.lon[7]):
                dsi.lon[7] = dsi.lon[6]
                print(f"{exp} {tid=}: changing lon=NaN at dtime=7 to {dsi.lon.data[6]:.3f}")
                ds[exp].loc[dict(id=tid)]['lon'] = dsi.lon
            if np.isnan(dsi.lon[8]):
                dsi.lon[8] = dsi.lon[9]
                print(f"{exp} {tid=}: changing lon=NaN at dtime=8 to {dsi.lon.data[9]:.3f}")
                ds[exp].loc[dict(id=tid)]['lon'] = dsi.lon

    return ds





def main():
    ds = open_track_dataset('../../tracker/data/')
    result = xr.DataArray([0,0,0], coords=[('x',[0,1,2])])
    print(f'{result=}')
    tasks = [
        add_data(result,1,1), 
        add_data(result,2,1),
        add_data(result,2,1)
    ]
    await asyncio.gather(*tasks)
    print(f'{result=}')

if __name__ == '__main__':
    asyncio.run(main())
