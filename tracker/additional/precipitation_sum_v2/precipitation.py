import sys
import multiprocessing as mp
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import xarray as xr
from src.analysis.load_SAIdata import open_mfdataset, Cases
import src.analysis.tracks as tracks


def read_CAM_prect(experiment, run_ensemble, stream_prec):
    case = Cases(f'hres.{experiment}.{run_ensemble}').select('atm',stream_prec)
    ds_cam = case.open_mfdataset(decode_times=False, cache=False)
    ds_cam = ds_cam.drop_vars(('date_written','time_written'))
    dt = (ds_cam.time[1] - ds_cam.time[0]).item() # time step (days)
    if 'time: mean' in getattr(ds_cam.PRECT, 'cell_methods', ''):
        ds_cam['time'] = ds_cam['time_bnds'].mean('nbnd', keep_attrs=False).assign_attrs(ds_cam.time.attrs)
    if dt == 0.125:
        print('converting 3hrly avg precipitation to 6hrly')
        ds_cam = ds_cam.coarsen(time=2).mean(keep_attrs=True)
    dt = (ds_cam.time[1] - ds_cam.time[0]).item() # time step (days)
    assert dt == 0.25, f'precipitation should be 6hrly, but {dt=}'
    return ds_cam


def get_timestamps(ds):
    """for each year, create array of timestamps at which TCs are active
    
    ds : xr.Dataset
        track dataset

    returns : 
    track_times : dict(str, arr)
        mapping from year to array of timestamps
    """
    track_times = {}
    for year in np.unique(ds.year):
        if np.isnan(year):
            continue
        dsi = ds.where(ds.year==year, drop=True)
        track_times_i = np.sort(np.unique(dsi.time))
        track_times_i = track_times_i[~np.isnan(track_times_i)]
        track_times[f'{int(year)}'] = track_times_i
        # print(int(year), len(track_times_i)/8) # days with active TCs
    return track_times


def circular_mask_TC(TC_lon, TC_lat, ds, r=5):
    """create circular mask of r degrees around TC
    using periodic boundary conditions for longitude

    TC_lon, TC_lat : float, float
        longitude and latitude of TC

    ds : xr.Dataset
        should contain lon and lat fields

    r : int | float
        radius of circular mask in degrees

    returns:
    circlemask : bool (ds.lat.size x ds.lon.size)
        true where distance to TC is less than or equal to r
    """
    assert r > 0, 'r must be positive'
    assert r < 30, 'can only handle radii <30degs to prevent issues near the poles'
    dlon = (ds.lon-TC_lon+180)%360 - 180 
    dlat = ds.lat-TC_lat
    circlemask = np.sqrt(dlat**2 + dlon**2) <= r
    return circlemask


def join_TC_masks(ids, ds_tracks, ds_cam, r=5):
    """create mask field that is true within r degrees of any TC

    ids : array of int (N x 2)
        point coordinates (id x dtime) of (N) TCs

    ds_tracks, ds_cam : xr.Dataset
        track / CAM dataset

    r : int | float
        radius of circular mask around individual TCs

    returns:
    joined_masks : bool (ds_cam.lat.size x ds_cam.lon.size)
        joined circular mask field
        true where distance to any TC is smaller than or equal to r
    """
    masks = [xr.zeros_like(ds_cam.lat * ds_cam.lon, dtype='bool')]
    for (tid,dt) in ids:  # loop through active TCs
        ds_TC = ds_tracks.isel(id=tid, dtime=dt)
        mask = circular_mask_TC(ds_TC.lon, ds_TC.lat, ds_cam, r)
        masks.append(mask)
    masks = xr.concat(masks, dim='n', coords='minimal', compat='override')
    joined_masks = masks.any('n')
    return joined_masks


def interpolate_time(time, da):
    """linearly interpolate da to time step 'time'

    time : float
        time to interpolate to

    da : xr.DataArray
        data to interpolate (must have time coordinate)

    returns:
    da_interp : xr.DataArray
        da interpolated to 'time'
    """
    if time in da.time:
        return da.sel(time=time)
    ax = da.dims.index('time')
    id1 = np.searchsorted(da.time, time)
    id0 = max(id1 - 1, 0)
    da_adj = da.isel(time=slice(id0,id1+1))
    if max(time-da_adj.time[-1], da_adj.time[0]-time) > 0.125:
        print(f"WARNING: extrapolating to {time=} from {da_adj.time.data}")
    interped = interp1d(da_adj.time, da_adj, axis=ax, fill_value='extrapolate')(time)
    da_interp = da.isel(time=0)
    da_interp.time.data = time
    da_interp.data = interped
    return da_interp


def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def main(experiment, run_ensemble):
    sys.stdout = open(f'log.{experiment}.{run_ensemble}', 'w', buffering=1)
    sys.stderr = sys.stdout
    exps = {'ref':'Reference', 'cnt':'RCP8.5', 'sai':'SAI2050'}
    Exp = exps[experiment]
    stream, stream_prec = ('h1','h2') if (experiment in ('ref','cnt')) and (run_ensemble <= 5) else ('h5','h3')
    savefile = f'PRECT_sum.{experiment}.{run_ensemble}.nc'
    print(f'{Exp=}, {stream=}, {stream_prec=}, {savefile=}')
    
    ds_cam = read_CAM_prect(experiment, run_ensemble, stream_prec)
    ds_tracks = tracks.load_tracks('../../jobs/Tracking_TC_RV.infext2/', ext='.infext')
    dsi_tracks = ds_tracks[Exp].where(ds_tracks[Exp].ens==run_ensemble, drop=True)  # select tracks for experiment and ensemble
    
    track_times = get_timestamps(dsi_tracks)
    precip_sum = {}
    dt = 3 * 60 * 60 # time difference (s) between track points
    for (year,timestamps) in track_times.items():  # loop through years
        print(year)
        pcip_sum_year = xr.zeros_like(ds_cam.lat*ds_cam.lon, dtype='float32')
        pcip_sum_year = (pcip_sum_year
                         .rename('PRECT')
                         .assign_coords(year=int(year))
                         .assign_attrs({'name':'PRECT','long_name':'precipitation sum','units':'mm/year'}))
        for time in timestamps[:]:  # loop through time stamps with active TCs
            if (time - timestamps[0])%1==0:
                print(f'current doy: {time-timestamps[0]+1}')
            PRECT_step = interpolate_time(time, ds_cam.PRECT) # 6hrly average precipitation (m/s)
            ids = np.nonzero(dsi_tracks.time == time).T  # IDs of active TCs at time step
            mask = join_TC_masks(ids, dsi_tracks, ds_cam, r=5) # mask (Nlat x Nlon), true within r degrees of any TC
            dPRECT = mask.astype('int') * PRECT_step * dt * 1000 # precipitation sum (mm) over 3hr interval between track points
            pcip_sum_year += dPRECT.compute()
        precip_sum[year] = pcip_sum_year.drop_vars(('time','id','dtime'))
    
    # save results
    ds = xr.concat(precip_sum.values(), dim='year').to_dataset()
    ds.to_netcdf(savefile)
    return 


if __name__ == '__main__':
    #for experiment in ['ref','cnt','sai']:
    #    for run_ensemble in range(1,7):
    #        mp.Process(target=main, args=(experiment, run_ensemble)).start()
    mp.Process(target=main, args=('ref',1)).start()
    mp.Process(target=main, args=('ref',5)).start()
    mp.Process(target=main, args=('cnt',3)).start()
    mp.Process(target=main, args=('cnt',5)).start()
    #main()
