"""Apply selection criteria and join TC tracks

Reads TC track data from 'directory' and CESM data from 'datadir',
Stitches TC tracks together that meet criteria
Results are stored in 'directory/outfile'

Call this script with
    >> python Tracking_TC_RV.py exp ens
where 
    exp : [ref, rcp, sai], high-res CESM experiment
    ens : [1,2,...,6], ensemble member
"""

import os
import shutil
import sys
import glob
import cftime
import numpy as np
import datetime
import math
import netCDF4 as netcdf
from scipy.interpolate import interp1d

experiment = sys.argv[1]  # run this script like python Tracking_TC_RV.py rcp 1
run_ensemble = int(sys.argv[2])
ROOT = '/home/jasperdj/TCs-under-SAI/tracker'
N_extend = 0 # number of steps to extend valid tracks using optional RV maxima (default 0, use -1 for unlimited)
SEEDS = False # set True to track TC seeds instead of TCs

if (experiment == 'ref') and (run_ensemble<=5):
    run_year_start  = 2002
    experiment_name = f'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_{run_year_start}-12.{run_ensemble:03d}'
    ddir            = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4/f02_t12/'+experiment_name+'/OUTPUT'
    datadir         = {'atm': ddir+'/atm/hist/3h', 
                       'atm_mon': ddir+'/atm/hist/monthly', 
                       'ocn':ddir+'/ocn/hist/daily'}
    stream          = 'h1' # file stream with 3 hourly instantaneous variables
    stream_prec     = 'h2' # file stream with 3 hourly average precipitation
    directory       = f'{ROOT}/data/RCP.started_{run_year_start}.{run_ensemble:03d}'
    gridfile        = f'{ROOT}/data/Atmosphere_0_25_DX_DY_AREA.nc'
    RVSEARCHRADIUS  = 200 # (km) search radius for TC at previous time step
    time_slice      = ('2003-01-01','2003-03-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.seeds.nc'
elif (experiment=='ref') and (run_ensemble==6):
    run_year_start  = 2002
    experiment_name = f'hres_b.e10.B2000_CAM5.f02_t12.started_{run_year_start}-12_without_SAI.001'
    ddir            = '/projects/0/nwo2021025/archive/'+experiment_name
    datadir         = {'atm': ddir+'/atm/hist',
                       'atm_mon': ddir+'/atm/hist',
                       'ocn':ddir+'/ocn/hist'}
    stream          = 'h5' # file stream with 3 hourly instantaneous variables
    stream_prec     = 'h3' # file stream with 6 hourly average precipitation
    directory       = f'{ROOT}/data/RCP.started_{run_year_start}.{run_ensemble:03d}'
    gridfile        = f'{ROOT}/data/Atmosphere_0_25_DX_DY_AREA.nc'
    RVSEARCHRADIUS  = 200 # (km) search radius for TC at previous time step
    time_slice      = ('2003-01-01','2008-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.seeds.nc'
elif (experiment == 'rcp') and (run_ensemble<=5):
    run_year_start  = 2092
    experiment_name = f'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_{run_year_start}-12.{run_ensemble:03d}'
    ddir            = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4/f02_t12/'+experiment_name+'/OUTPUT'
    datadir         = {'atm': ddir+'/atm/hist/3h', 
                       'atm_mon': ddir+'/atm/hist/monthly', 
                       'ocn':ddir+'/ocn/hist/daily'}
    stream          = 'h1' # file stream with 3 hourly instantaneous variables
    stream_prec     = 'h2' # file stream with 3 hourly average precipitation
    directory       = f'{ROOT}/data/RCP.started_{run_year_start}.{run_ensemble:03d}'
    gridfile        = f'{ROOT}/data/Atmosphere_0_25_DX_DY_AREA.nc'
    RVSEARCHRADIUS  = 200 # (km) search radius for TC at previous time step
    time_slice      = ('2093-01-01','2098-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.seeds.nc'
elif (experiment=='rcp') and (run_ensemble==6):
    run_year_start  = 2092
    experiment_name = f'hres_b.e10.B2000_CAM5.f02_t12.started_{run_year_start}-12_without_SAI.001'
    ddir            = '/projects/0/nwo2021025/archive/'+experiment_name
    datadir         = {'atm': ddir+'/atm/hist',
                       'atm_mon': ddir+'/atm/hist',
                       'ocn':ddir+'/ocn/hist'}
    stream          = 'h5' # file stream with 3 hourly instantaneous variables
    stream_prec     = 'h3' # file stream with 6 hourly average precipitation
    directory       = f'{ROOT}/data/RCP.started_{run_year_start}.{run_ensemble:03d}'
    gridfile        = f'{ROOT}/data/Atmosphere_0_25_DX_DY_AREA.nc'
    RVSEARCHRADIUS  = 200 # (km) search radius for TC at previous time step
    time_slice      = ('2093-01-01','2098-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.seeds.nc'
elif experiment == 'sai':
    run_year_start  = 2092
    experiment_name = 'hres_b.e10.B2000_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)
    ddir            = '/projects/0/nwo2021025/archive/'+experiment_name
    datadir         = {'atm': ddir+'/atm/hist', 
                       'atm_mon': ddir+'/atm/hist', 
                       'ocn':ddir+'/ocn/hist'}
    stream          = 'h5' # file stream with 3 hourly instantaneous variables
    stream_prec     = 'h3' # file stream with 6 hourly average precipitation
    directory       = f'{ROOT}/data/SAI.started_'+str(run_year_start)+'.'+str(run_ensemble).zfill(3)
    gridfile        = f'{ROOT}/data/Atmosphere_0_25_DX_DY_AREA.nc'
    RVSEARCHRADIUS  = 200 # (km) search radius for TC at previous time step
    time_slice      = ('2093-01-01','2098-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.seeds.nc'
# elif experiment == 'mres': # SAI mres data
#   experiment_name = 'mres_b.e10.B2000_CAM5.f05_t12.001'
#   ddir            = '/projects/0/nwo2021025/archive/'+experiment_name
#   datadir         = {'atm': ddir+'/atm/hist/',
#                       'atm_mon': ddir+'/atm/hist/',
#                       'ocn':ddir+'/ocn/hist/'}
#   stream          = 'h4' # 6-hourly averages
#   stream_prec     = 'h4' # 6-hourly average precipitation
#   stream3D        = 'h3' # only used for mres to calculate U250, V250 and T850
#   directory       = '{ROOT}/data/SAI.mres'
#   gridfile        = '{ROOT}/data/Atmosphere_0_5_DX_DY_AREA.nc'
#   RVSEARCHRADIUS = 400 # (km) search radius for TC at previous time step
#   NOUT            = 4 # number of time steps per day
#   outfile         = 'TC_tracker_results.{experiment}.nc'
else:
    raise ValueError(f'Incorrect options provided ({experiment=}, {run_ensemble=}) provided, experiment should be in [ref,rcp,sai], ensemble in [1-6].')

outdir = directory  # directory for output file
files = sorted(glob.glob(f'{datadir["atm"]}/{experiment_name}.cam2.{stream}.*.nc')) # input (CAM) files
# files = files[6:373]
outdir = '.'
outfile = 'test.reg.2mon.nc'
files = files[6:20]

# def chunk_data(files):
#   """divide data into daily chunks"""
#   n = 0 
#   dt = 1/NOUT # output time step in days
#   with netcdf.Dataset(files[0],'r') as fh:
#       time = fh.variables['time'][:]
#       if 'time: mean' in getattr(fh['U850'], 'cell_methods', ''):
#           time = time - 0.5*dt
#       dates = [cftime.num2date(t, fh['time'].units, fh['time'].calendar).strftime("%Y%m%d") for t in time]
#       chunks = [[dates[0],[0],0,0]]
#   for fid in range(len(files)):
#       with netcdf.Dataset(files[fid],'r') as fh:
#           time = fh.variables['time'][:]
#           if 'time: mean' in getattr(fh['U850'], 'cell_methods', ''):
#               time = time - 0.5*dt
#           dates = [cftime.num2date(t, fh['time'].units, fh['time'].calendar).strftime("%Y%m%d") for t in time]
#           for i,newdate in enumerate(dates):
#               (olddate, fids0, i0, n0) = chunks[-1]
#               if newdate != olddate:
#                   if fids0[-1] != fid and i!=0:
#                       chunks[-1][1].append(fid)
#                   chunks[-1][3] = i0 + n - n0
#                   chunks.append([newdate,[fid],i,n])
#               n += 1
#   (olddate, fids0, i0, n0) = chunks[-1]
#   chunks[-1][3] = i0 + n - n0
#   return chunks


def chunk_data(files):
    """divide data into daily chunks
    note: date is determined by center of time_bnds, not time!
    
    result: chunks = [(date, fileids, t0, t1), ...]
        such that netcdf.MFDataset(files[fileids])['time'][t0:t1] gives
        the timestamps corresponding to date
    """
    print("start chunking")
    n = 0 # counter for time step
    with netcdf.Dataset(files[0],'r') as fh:
        time = fh['time_bnds'][:].mean(axis=-1)
        dates = [cftime.num2date(t, fh['time'].units, fh['time'].calendar).strftime("%Y-%m-%d") for t in time]
        chunks = [[dates[0],[0],0,0]]
    for fid in range(len(files)):
        with netcdf.Dataset(files[fid],'r') as fh:
            time = fh['time_bnds'][:].mean(axis=-1)
            dates = [cftime.num2date(t, fh['time'].units, fh['time'].calendar).strftime("%Y-%m-%d") for t in time]
            for i,newdate in enumerate(dates):
                (olddate, fids0, i0, n0) = chunks[-1]
                if newdate != olddate:
                    if fids0[-1] != fid and i!=0:
                        chunks[-1][1].append(fid)
                    chunks[-1][3] = i0 + n - n0
                    chunks.append([newdate,[fid],i,n])
                n += 1
    (olddate, fids0, i0, n0) = chunks[-1]
    chunks[-1][3] = i0 + n - n0
    print(f"created {len(chunks)} daily time chunks")
    return chunks
    
    
def getattstr(dataset, name):
    ncvar = dataset[name]
    attrs = {a:getattr(ncvar,a) for a in ncvar.ncattrs() if a not in ['_FillValue','missing_value']}
    if name == 'PSL':
        attrs['units'] = 'h'+attrs['units']
        if 'cell_methods' in attrs and 'time: mean' in attrs['cell_methods']:
            name = 'PSLmon' 
    if name == 'SST':
        attrs['missing_value'] = 100
    return f"{name}:{attrs}"


def ReadinData(filenames, x_diff, y_diff, t1, t2):
    tmpfiles = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(tmpdir,'*.nc')))]
    for file in filenames:
        if os.path.basename(file) not in tmpfiles:  # copy nc file to tmpdir if not present
            print(f"copying {os.path.basename(file)} to tmpdir")
            shutil.copyfile(file, os.path.join(tmpdir,os.path.basename(file)))
    for file in tmpfiles:
        if (file not in [os.path.basename(f) for f in filenames]) and (f'.{stream}.' in file):
            print(f"remove {file} from tmpdir")
            os.remove(os.path.join(tmpdir,file)) # delete unused nc files from tmpdir
    tmpfilenames = [os.path.join(tmpdir, os.path.basename(f)) for f in filenames]
    print(f"reading {[os.path.basename(f) for f in tmpfilenames]} in tmpdir, time steps {t1}-{t2}")
    with netcdf.MFDataset(tmpfilenames, 'r', aggdim='time') as fh:
        time        = fh.variables['time'][t1:t2]   # Time (end of interval)
        if 'time: mean' in getattr(fh['U850'], 'cell_methods', ''):
            time = fh['time_bnds'][t1:t2].mean(axis=-1) # Center time if data represents time average
        timestamp   = cftime.num2date(time, fh['time'].units, fh['time'].calendar)
        date        = fh['time_bnds'][t1:t2].mean(axis=-1) # always center time for date
        date        = np.array([cftime.num2date(d, fh['time'].units, fh['time'].calendar).strftime('%Y-%m-%d') for d in date])
        lon         = fh.variables['lon'][:]                       #Longitude
        lat         = fh.variables['lat'][126:642]      #Latitude
        lat_weight  = fh.variables['gw'][127:641]
        pres        = fh.variables['PSL'][t1:t2, 127:641] / 100.0   #Sea level pressure (hPa)
        U_10        = fh.variables['U10'][t1:t2, 127:641]   #10-meter wind speed (m/s)
        u_vel_850   = fh.variables['U850'][t1:t2, 126:642]  #Zonal velocity at 850 hPa (m/s)
        v_vel_850   = fh.variables['V850'][t1:t2, 126:642]  #Meridional velocity at 850 hPa (m/s)
        u_vel_250   = fh.variables['U250'][t1:t2, 127:641]  #Zonal velocity at 250 hPa (m/s)
        v_vel_250   = fh.variables['V250'][t1:t2, 127:641]  #Meridional velocity at 250 hPa (m/s)
        if stream_prec == stream:
            prec    = fh.variables['PRECT'][t1:t2, 127:641]
            if track_attrs[11] == '':
                track_attrs[11] = getattstr(fh, 'PRECT')
        else:
            prec    = ReadinDataPRECT(timestamp)
        
        if track_attrs[0] == '':
            track_attrs[0] = getattstr(fh, 'time')
            track_attrs[1] = getattstr(fh, 'lon')
            track_attrs[2] = getattstr(fh, 'lat')
            track_attrs[3] = getattstr(fh, 'PSL')
            track_attrs[5] = getattstr(fh, 'U10')
            track_attrs[9]  = "RV:{'units':'1/s', 'long_name':'relative vorticity'}"
            track_attrs[10] = f"Vshear:{{'units':'{fh['U850'].units}', 'long_name':'vertical wind shear 250-850hPa'}}"
            

    #Add periodic boundaries
    lon_2, u_vel_850    = PeriodicBoundaries3D(time, lon, lat, u_vel_850)
    lon, v_vel_850      = PeriodicBoundaries3D(time, lon, lat, v_vel_850)

    #Determine the vorticity at 850 and 250 hPa
    vor_850         = np.ma.masked_all((len(time), len(lat) - 2, len(lon) - 2))

    for time_i in range(len(time)):
        #Vorticity for the two level
        vor_850[time_i] = RelativeVorticity(u_vel_850[time_i], v_vel_850[time_i], x_diff, y_diff)
    
    #Get the correct lon/lat dimensions
    lon, lat        = lon[1:-1], lat[1:-1]

    return time, timestamp, date, lon, lat, lat_weight, pres, U_10, vor_850, u_vel_850[:, 1:-1, 1:-1], v_vel_850[:, 1:-1, 1:-1], u_vel_250, v_vel_250, prec

    
# # original function
#
# def ReadinDataPRECT(timestamp):
#   """Get precipitation from separate file stream, linearly interpolating to timestamp"""
#   files = sorted(glob.glob(f'{datadir["atm"]}/{experiment_name}.cam2.{stream_prec}.*.nc'))
#   dates = [file[-19:-3] for file in files]
#   date0, date1 = timestamp[0].strftime('%Y-%m-%d'), timestamp[-1].strftime('%Y-%m-%d')
#   fid0 = max(np.searchsorted(dates, date0)-1,0)
#   fid1 = np.searchsorted(dates, date1)+1
#   with netcdf.MFDataset(files[fid0:fid1], 'r') as fh: 
#       nc_time = fh['time'][:] # Time (end of interval)
#       if 'time: mean' in getattr(fh['PRECT'], 'cell_methods', ''):
#           nc_time = fh['time_bnds'][:].mean(axis=-1) # Center time if data represents time average
#       times = cftime.date2num(timestamp, fh['time'].units, fh['time'].calendar) # convert timestamp to same units
#       tids = [np.argmin(np.abs(nc_time-t)) for t in times]
#       prec = fh['PRECT'][tids,127:641]
#       if times[0] < nc_time[0] or times[-1] > nc_time[-1]:
#           ids_invalid = (times < nc_time[0]) | (times > nc_time[-1])
#           prec[ids_invalid] = np.nan
#           print(f"WARNING: Failed to interpolate PRECT to {times[ids_invalid]=}, available times: {nc_time=}")
#       if track_attrs[11] == '':
#           track_attrs[11] = getattstr(fh, 'PRECT')
#   return prec


#def ReadinDataPRECT(timestamp):
#    """Get precipitation from separate file stream, linearly interpolating to timestamp"""
#    files = sorted(glob.glob(f'{datadir["atm"]}/{experiment_name}.cam2.{stream_prec}.*.nc'))
#    dates = [file[-19:-3] for file in files]
#    date0, date1 = timestamp[0].strftime('%Y-%m-%d'), timestamp[-1].strftime('%Y-%m-%d')
#
#    # first obtain three files that should contain all steps in timestamp
#    fid0 = max(np.searchsorted(dates, date0)-1,0)
#    fid1 = np.searchsorted(dates, date1)+1
#    tmpfiles = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(tmpdir,'*.nc')))]
#    for file in files[fid0:fid1]:
#        if os.path.basename(file) not in tmpfiles: # copy nc file to tmpdir if not present
#            print(f"copying {os.path.basename(file)} to tmpdir")
#            shutil.copyfile(file, os.path.join(tmpdir,os.path.basename(file)))
#    for file in tmpfiles:
#        if (file not in [os.path.basename(f) for f in files[fid0:fid1]]) and (f'.{stream_prec}.' in file):
#            print(f"remove {file} from tmpdir")
#            os.remove(os.path.join(tmpdir,file)) # delete unused nc files from tmpdir
#    tmpfilenames = [os.path.join(tmpdir, os.path.basename(f)) for f in files[fid0:fid1]]
#
#    # loading all data from the three files is slow, so first only check time and remove the redundant files
#    for file in tmpfilenames.copy():
#        with netcdf.Dataset(file, 'r') as fh:
#            times = cftime.date2num(timestamp, fh['time'].units, fh['time'].calendar) # convert timestamp to same units
#            t0,t1 = fh['time'][0], fh['time'][-1]
#        if t1 < times[0]:
#            tmpfilenames.remove(file)
#            continue
#        if t0 > times[-1]:
#            tmpfilenames.remove(file)
#            continue
#
#    # load the data
#    print(f"reading PRECT from {tmpfilenames}")
#    with netcdf.MFDataset(tmpfilenames, 'r') as fh: 
#        nc_time = fh['time'][:] # Time (end of interval)
#        nc_tbnd = fh['time_bnds'][:] # time bounds
#        nc_prec = fh['PRECT'][:,127:641] # precipitation
#        if 'time: mean' in getattr(fh['PRECT'], 'cell_methods', ''):
#            nc_time = nc_tbnd.mean(axis=-1) # Center time if data represents time average
#            if np.diff(nc_tbnd[0]) == 0.125: # average 3hrly data to 6hrly
#                nc_time = np.mean([nc_time[1::2],nc_time[:-1:2]], axis=0)
#                nc_prec = np.mean([nc_prec[1::2],nc_prec[:-1:2]], axis=0)
#        times = cftime.date2num(timestamp, fh['time'].units, fh['time'].calendar) # convert timestamp to same units
#        print(f"interpolating PRECT from {len(nc_prec)=} to {len(nc_time)=} timestamps")
#        prec = interp1d(nc_time, nc_prec, axis=0, bounds_error=False)(times)
#        if times[0] < nc_time[0] or times[-1] > nc_time[-1]:
#            ids_invalid = (times < nc_time[0]) | (times > nc_time[-1])
#            prec[ids_invalid] = np.nan
#            print(f"\nWARNING: Failed to interpolate PRECT to {times[ids_invalid]=} "
#                + f"(={[str(t) for t in timestamp[ids_invalid]]}), available times: {nc_time=}")
#        if track_attrs[11] == '':
#            track_attrs[11] = getattstr(fh, 'PRECT')
#    return prec
    

def ReadinDataPRECT(timestamp):
    """Get precipitation from separate file stream, linearly interpolating to timestamp"""
    files = sorted(glob.glob(f'{datadir["atm"]}/{experiment_name}.cam2.{stream_prec}.*.nc'))
    dates = [file[-19:-3] for file in files]
    date0, date1 = timestamp[0].strftime('%Y-%m-%d'), timestamp[-1].strftime('%Y-%m-%d')

    # first obtain three files that should contain all steps in timestamp
    fid0 = max(np.searchsorted(dates, date0)-1,0)
    fid1 = np.searchsorted(dates, date1)+1
    tmpfiles = [os.path.basename(f) for f in sorted(glob.glob(os.path.join(tmpdir,'*.nc')))]
    for file in files[fid0:fid1]:
        if os.path.basename(file) not in tmpfiles: # copy nc file to tmpdir if not present
            print(f"copying {os.path.basename(file)} to tmpdir")
            shutil.copyfile(file, os.path.join(tmpdir,os.path.basename(file)))
    for file in tmpfiles:
        if (file not in [os.path.basename(f) for f in files[fid0:fid1]]) and (f'.{stream_prec}.' in file):
            print(f"remove {file} from tmpdir")
            os.remove(os.path.join(tmpdir,file)) # delete unused nc files from tmpdir
    tmpfilenames = [os.path.join(tmpdir, os.path.basename(f)) for f in files[fid0:fid1]]

    # loading all data from the three files is slow, so first only check time and remove the redundant files
    for file in tmpfilenames.copy():
        with netcdf.Dataset(file, 'r') as fh:
            times = cftime.date2num(timestamp, fh['time'].units, fh['time'].calendar) # convert timestamp to same units
            t0,t1 = fh['time'][0], fh['time'][-1]
        if t1 < times[0]:
            tmpfilenames.remove(file)
            continue
        if t0 > times[-1]:
            tmpfilenames.remove(file)
            continue

    # load the data
    print(f"reading PRECT from {[os.path.basename(f) for f in tmpfilenames]} in tmpdir")
    with netcdf.MFDataset(tmpfilenames, 'r') as fh:
        nc_time = fh['time'][:] # Time (end of interval)
        nc_tbnd = fh['time_bnds'][:] # time bounds
        nc_prec = fh['PRECT'][:,127:641] # precipitation
        if 'time: mean' in getattr(fh['PRECT'], 'cell_methods', ''):
            nc_ctime = nc_tbnd.mean(axis=-1) # Center time if data represents time average
            if np.diff(nc_tbnd[0]) == 0.125: # average 3hrly data to 6hrly
                nc_ctime = np.mean([nc_ctime[1::2],nc_ctime[:-1:2]], axis=0)
                print(f"averaging 3hrly PRECT ({len(nc_time)} steps) to 6hrly ({len(nc_ctime)} steps)")
                nc_prec = np.mean([nc_prec[1::2],nc_prec[:-1:2]], axis=0)
        times = cftime.date2num(timestamp, fh['time'].units, fh['time'].calendar) # convert timestamp to same units
        print(f"interpolating PRECT from {len(nc_ctime)} to {len(times)} timestamps")
        prec = interp1d(nc_ctime, nc_prec, axis=0, bounds_error=False)(times)
        if times[0] < nc_tbnd[0,0] or times[-1] > nc_tbnd[-1,1]:
            ids_invalid = (times < nc_tbnd[0,0]) | (times > nc_tbnd[-1,1])
            prec[ids_invalid] = np.nan
            print(f"\nWARNING: Failed to interpolate PRECT to {times[ids_invalid]=} "
                + f"(={[str(t) for t in timestamp[ids_invalid]]}), available times: {nc_time=}")
        if track_attrs[11] == '':
            track_attrs[11] = getattstr(fh, 'PRECT')
    return prec


def ReadinDataPSLMonth(month, year):
    """Read in the sea-level pressure (monthly)"""
    filename = f'{datadir["atm_mon"]}/{experiment_name}.cam2.h0.{year}-{month:02}.nc'
    print(f"Reading {filename}")
    fh = netcdf.Dataset(filename, 'r')
    pres_month      = fh.variables['PSL'][0, 127:641] / 100.0   #Sea level pressure (hPa)
    if track_attrs[4] == '':
        track_attrs[4] = getattstr(fh, 'PSL')
    fh.close()
    return pres_month


def ReadinDataSST(month, year):
    """Read in the Sea surface temperature (daily averages)"""
    filename = f'{datadir["ocn"]}/{experiment_name}.pop.h.nday1.{year}-{month:02}-01.nc'
    try:
        fh = netcdf.Dataset(filename, 'r')
        print(f"Reading {filename}")
        if 'time: mean' in getattr(fh['SST'], 'cell_methods', ''):
            time_ocn = fh['time_bound'][:].mean(axis=-1)
        else:
            time_ocn = fh['time'][:]
        time_ocn    = cftime.num2date(time_ocn, fh['time'].units, fh['time'].calendar)
        lon_ocn     = fh.variables['TLONG'][:]
        lat_ocn     = fh.variables['TLAT'][:]
        KMT     = fh.variables['KMT'][:]
        SST         = fh.variables['SST'][:]    #Sea surface temperature
        if track_attrs[7] == '':
            track_attrs[7] = getattstr(fh, 'SST')
            track_attrs[8] = f"SSTmon:{{'units':'{fh['SST'].units}'}}"
        fh.close()
        for time_i in range(len(time_ocn)): #Mask all land elements
            SST[time_i] = np.ma.masked_where(KMT == -1.0, SST[time_i])
            lon_ocn = np.ma.masked_where(KMT == -1.0, lon_ocn)
            lat_ocn = np.ma.masked_where(KMT == -1.0, lat_ocn)
    except FileNotFoundError:
        print(f"WARNING: Could not find ocean file {filename}")
        if track_attrs[7] == '':
            track_attrs[7] = "SST:{'units':'degC'}"
            track_attrs[8] = "SSTmon:{'units':'degC'}"
        return None, None, None, None, None
    return time_ocn, lon_ocn, lat_ocn, SST, np.mean(SST, axis = 0)


def ReadinData_RV_Max(RVfile, timestamp):
    """Read in the PSL minima for tracking"""
    try:
        HEAT_data = netcdf.Dataset(RVfile, 'r')
    except FileNotFoundError:
        print(f"Could not find {RVfile}.")
        raise

    time_RV     = HEAT_data.variables['time']
    time_RV     = cftime.num2date(time_RV[:], time_RV.units, time_RV.calendar)
    try:
        time_index  = np.where(timestamp == time_RV)[0][0]
    except IndexError:
        print(f"ERROR: Could not find {timestamp=} in {RVfile} where {time_RV=}")
        raise
    lon_index   = HEAT_data.variables['lon_index'][time_index] 
    lat_index   = HEAT_data.variables['lat_index'][time_index]  
    temp_anom   = HEAT_data.variables['TEMP'][time_index]
    if track_attrs[12] == '':
        track_attrs[12] = getattstr(HEAT_data, 'TEMP').replace('longname','long_name')

    HEAT_data.close()

    #Only retain the non-masked elements
    mask_index  = np.where(lon_index.mask == True)[0][0]
    lon_index   = lon_index[:mask_index]
    lat_index   = lat_index[:mask_index]
    temp_anom   = temp_anom[:mask_index]

    return lon_index.astype(int), lat_index.astype(int), temp_anom


def PeriodicBoundaries2D(lon, lat, field, lon_grids = 1):
    """Add periodic zonal boundaries for 2D field"""

    #Empty field with additional zonal boundaries
    lon_2           = np.zeros(len(lon) + lon_grids * 2)
    field_2         = np.ma.masked_all((len(lat), len(lon_2)))
    
    #Get the left boundary, which is the right boundary of the original field
    lon_2[:lon_grids]   = lon[-lon_grids:] - 360.0
    field_2[:, :lon_grids]  = field[:, -lon_grids:]

    #Same for the right boundary
    lon_2[-lon_grids:]  = lon[:lon_grids] + 360.0
    field_2[:, -lon_grids:] = field[:, :lon_grids]

    #And the complete field
    lon_2[lon_grids:-lon_grids]     = lon
    field_2[:, lon_grids:-lon_grids]    = field

    return lon_2, field_2   


def PeriodicBoundaries3D(time, lon, lat, field, lon_grids = 1):
    """Add periodic zonal boundaries for 3D field"""

    #Empty field with additional zonal boundaries
    lon_2               = np.zeros(len(lon) + lon_grids * 2)
    field_2             = np.ma.masked_all((len(time), len(lat), len(lon_2)))
    
    #Get the left boundary, which is the right boundary of the original field
    lon_2[:lon_grids]       = lon[-lon_grids:] - 360.0
    field_2[:, :, :lon_grids]   = field[:, :, -lon_grids:]

    #Same for the right boundary
    lon_2[-lon_grids:]      = lon[:lon_grids] + 360.0
    field_2[:, :, -lon_grids:]  = field[:, :, :lon_grids]

    #And the complete field
    lon_2[lon_grids:-lon_grids]     = lon
    field_2[:, :, lon_grids:-lon_grids]     = field

    return lon_2, field_2   


def RelativeVorticity(u_vel, v_vel, x_diff, y_diff):
    """Determines the relative vorticity of the field"""

    #Take the meridional difference of the zonal wind
    u_vel_diff  = u_vel[2:] - u_vel[:-2]
    u_vel_diff  = u_vel_diff[:, 1:-1]

    #Take the zonal difference of the meridional wind
    v_vel_diff  = v_vel[:, 2:] - v_vel[:, :-2]
    v_vel_diff  = v_vel_diff[1:-1]

    #Determine the relative vorticity
    vorticity   = (v_vel_diff / x_diff) - (u_vel_diff / y_diff)

    return vorticity


def FieldPartition(lon, lat, lon_index, lat_index, zonal_extent, meridional_extent, lat_weight, field):
    """Gets a part of the field and reshapes near 0E or 360E boundary"""
    
    #Get the field around given coordinate
    lon_west    = lon_index - zonal_extent
    lon_east    = lon_index + zonal_extent + 1
    lat_south   = lat_index - meridional_extent
    lat_north   = lat_index + meridional_extent + 1

    if lat_south < 0:
        #Southern boundary is reached, set index to zero
        lat_south = 0

    #Get the dimensions for the lonxlat field
    lat_field   = lat[lat_south:lat_north]
    weight      = lat_weight[lat_south:lat_north]
    weight      = weight / np.sum(weight)
    lon_field   = np.zeros(lon_east - lon_west)
    new_field   = np.ma.masked_all((len(lat_field), len(lon_field)))

    if lon_east > len(lon):
        #Eastern part is on eastern hemisphere
        lon_1               = lon[lon_west:]
        lon_2               = lon[:lon_east - len(lon)] + 360.0
        lon_field[:len(lon_1)]      = lon_1
        lon_field[len(lon_1):]      = lon_2
        new_field[:, :len(lon_1)]   = field[lat_south:lat_north, lon_west:]
        new_field[:, len(lon_1):]   = field[lat_south:lat_north, :lon_east - len(lon)]

    elif lon_west < 0:
        #Western part is on western hemisphere
        lon_1               = lon[lon_west:] - 360.0
        lon_2               = lon[:lon_east]
        lon_field[:len(lon_1)]      = lon_1
        lon_field[len(lon_1):]      = lon_2
        new_field[:, :len(lon_1)]   = field[lat_south:lat_north, lon_west:]
        new_field[:, len(lon_1):]   = field[lat_south:lat_north, :lon_east]

    else:
        #Normal field can be retained
        lon_field   = lon[lon_west:lon_east]
        new_field   = field[lat_south:lat_north, lon_west:lon_east]

    return lon_field, lat_field, weight, new_field


def RadiusMaximumWind(lon, lat, lon_index, lat_index, vel_field_all, reduced_radius = False):
    """Determines the radius of maximum wind speeds near pressure minimum (eye of TC)
    Radius is determined in four quadrants and averaged to smooth out fixed grid distances
    Returns radius of maximum wind speeds (km) and the maximum wind speeds (m/s)"""

    #Get 4x4 grid around given coordinate
    lon_field, lat_field, weight, vel_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, vel_field_all)

    if reduced_radius:
        #Radius is too large, reduce region to 2 x 2
        lon_field, lat_field, weight, vel_field = FieldPartition(lon, lat, lon_index, lat_index, 7, 9, lat_weight, vel_field_all)
    
    # compute radius of max. velocity in four quadrants
    (y0,x0) = (shp//2 for shp in vel_field.shape)
    distances = [np.nan] * 4
    for q,quadrant in enumerate((
        (slice(0,y0),slice(0,x0+1)),   # NW
        (slice(y0,None),slice(0,x0)),  # SW
        (slice(y0+1,None),slice(x0,None)), # SE
        (slice(0,y0+1),slice(x0+1,None)),  # NE
    )):
        lon_field_q = lon_field[quadrant[1]]
        lat_field_q = lat_field[quadrant[0]]
        vel_field_q = vel_field[quadrant]
        
        #Find the minimum index
        index_max   = np.where(vel_field_q == np.max(vel_field_q))

        #Get the index of the local minimum of current field 
        lon_vel_max = lon_field_q[index_max[1][0]]
        lat_vel_max = lat_field_q[index_max[0][0]]

        #Determine distance (km) between pressure minima and velocity maxima
        distances[q]    = Distance(lon[lon_index], lat[lat_index], lon_vel_max, lat_vel_max) / 1000.0
    
    distance = np.nanmean(distances)
    
    if track_attrs[6] == '':
        track_attrs[6]  = ("RMW:{'units':'km', 'long_name':'radius of maximum wind'," 
            + "'description':'average RMW of four quadrants in 4x4deg area (2x2 if RMW would be more than 400km)'}")

    return np.max(vel_field), distance


def VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850, v_vel_850, u_vel_250, v_vel_250):
    """Determines the vertical wind shear near pressure minimum (eye of TC)"""

    #Get 8x8 grid, centered at given coordinate
    lon_field, lat_field, weight, u_vel_850_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, u_vel_850)
    lon_field, lat_field, weight, v_vel_850_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, v_vel_850)
    lon_field, lat_field, weight, u_vel_250_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, u_vel_250)
    lon_field, lat_field, weight, v_vel_250_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, v_vel_250)

    #Get weighted field
    area    = np.zeros(u_vel_850_field.shape)

    for lat_i in range(len(lat_field)):
        area[lat_i] = weight[lat_i]
    
    #Get 3 degrees from center  
    x, y    = np.meshgrid(lon_field, lat_field)
    rad_deg = np.sqrt((x - lon[lon_index])**2.0 + (y - lat[lat_index])**2.0)
    area    = np.ma.masked_where(np.sqrt((x - lon[lon_index])**2.0 + (y - lat[lat_index])**2.0) > 3.0, area)

    #Normalise to total
    area    = area / np.sum(area)

    #Take the spatial average
    u_vel_850_field = np.sum(u_vel_850_field * area)
    v_vel_850_field = np.sum(v_vel_850_field * area)
    u_vel_250_field = np.sum(u_vel_250_field * area)
    v_vel_250_field = np.sum(v_vel_250_field * area)
    
    u_vel_shear = np.abs(u_vel_250_field - u_vel_850_field)
    v_vel_shear = np.abs(v_vel_250_field - v_vel_850_field)
    vel_shear   = np.sqrt(u_vel_shear**2.0 + v_vel_shear**2.0)

    return vel_shear


def MaximumPrecipitation(lon, lat, lon_index, lat_index, prec_field):
    """Determines the maximum precipitation near pressure minimum (eye of TC)
    Returns maximum precipitation rate (mm / day)"""

    #Get 4x4 grid around given coordinate
    lon_field, lat_field, weight, prec_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, prec_field)

    return np.max(prec_field)


def SeaSurfaceTemperatureTC(lon_ocn, lat_ocn, SST_day, SST_month, lon_pres, lat_pres):
    """Returns the SST of the pressure minima in the TC"""
    
    if SST_day is None or SST_month is None:
        return 100, 100 # if no SST available, fill with 100 (fulfills SST>25 requirement)

    #Retain the corresponding lat for the ocean
    index       = np.where(np.abs(lat_ocn - lat_pres) < 0.2)
    index_lat   = index[0]
    index_lon   = index[1]
    lon     = lon_ocn[index]

    #Retain the corresponding lon for the ocean
    index       = np.where(np.abs(lon - lon_pres) < 0.3)
    index_lat   = index_lat[index]
    index_lon   = index_lon[index]

    if len(index_lat) == 0:
        #Somewhere above masked land
        return SST_day[0, 0], SST_month[0, 0]

    #Empty array for determining the distance
    distance    = np.zeros(len(index_lat)) + 1000.0

    for index_i in range(len(index_lat)):
        #Get the corresponding coordinates
        lon = lon_ocn[index_lat[index_i], index_lon[index_i]]
        lat = lat_ocn[index_lat[index_i], index_lon[index_i]]

        #Distance in km
        distance[index_i]   = Distance(lon, lat, lon_pres, lat_pres) / 1000.0

    #Determine the minimum distance
    index_min   = np.argmin(distance)
    index_lat   = index_lat[index_min]
    index_lon   = index_lon[index_min]

    SST_TC_day  = SST_day[index_lat, index_lon]
    SST_TC_month    = SST_month[index_lat, index_lon]
    
    return SST_TC_day, SST_TC_month


def PressureRVMaxima(lon, lat, lon_index, lat_index, vor_field):
    """Determines the RV Maxima near the sea-level pressure minima"""

    #Get 4x4 grid around given coordinate (add -1 for SH)
    lon_field, lat_field, weight, vor_field = FieldPartition(lon, lat, lon_index, lat_index, 3, 4, lat_weight, vor_field)
    vor_field               = vor_field * np.sign(lat[lat_index])

    #Return maximum vorticity
    return np.max(vor_field)


def Distance(lon_1, lat_1, lon_2, lat_2):
    """Returns distance (m) of two points located at the globe coordinates need input in degrees"""

    #Convert to radians
    lon_1, lat_1, lon_2, lat_2 = np.radians([lon_1, lat_1, lon_2, lat_2]) 

    #Haversine formula 
    d_lon   = lon_2 - lon_1 
    d_lat   = lat_2 - lat_1 
    a   = math.sin(d_lat/2.0)**2 + math.cos(lat_1) * math.cos(lat_2) * math.sin(d_lon/2.0)**2
    c   = 2.0 * math.asin(np.sqrt(a)) 
    r   = 6371000.0 # Radius of earth in meters
    
    return c * r #Distance between two points in meter


def out_of_latbounds(lat_index, lat):
    """Returns (bool) whether a TC candidate has reached 60N/S latitude"""
    return (lat_index == 0) or (lat_index == len(lat) - 1)


def has_no_closed_circulation(lon, lat, lon_index, lat_index, lat_weight, u_vel_850_i, v_vel_850_i):
    """Check if the wind at 850 hPa is cyclonic and exceeds the minimum wind requirement (7.5 m/s) in four quadrants within a box of 8x8 degrees around the RV maximum"""
    #Get a part of the field
    lon_field, lat_field, weight, u_vel_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, u_vel_850_i)
    lon_field, lat_field, weight, v_vel_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, v_vel_850_i)

    #Determine the spatial average over the field
    u_vel_mean  = np.mean(u_vel_field, axis = 1)
    v_vel_mean  = np.mean(v_vel_field, axis = 1)
    u_vel_mean  = np.sum(u_vel_mean * weight)
    v_vel_mean  = np.sum(v_vel_mean * weight)

    #Remove the time mean to remove background flow
    u_vel_field = u_vel_field - u_vel_mean
    v_vel_field = v_vel_field - v_vel_mean

    #Get the indices of pressure minimum
    lon_index   = (np.abs(lon_field - lon[lon_index])).argmin()
    lat_index   = (np.abs(lat_field - lat[lat_index])).argmin()

    #Take section through eye of cyclone
    lon_index_east  = lon_index + 3
    lon_index_west  = lon_index - 3
    lat_index_north = lat_index + 4
    lat_index_south = lat_index - 4

    if lat_index_south < 0:
        lat_index_south = 0

    #Retain the velocity (times -1 for SH)
    u_vel_north = u_vel_field[lat_index+1:lat_index_north+1, lon_index] * np.sign(lat_field[lat_index])
    u_vel_south = u_vel_field[lat_index_south:lat_index, lon_index] * np.sign(lat_field[lat_index])
    v_vel_west  = v_vel_field[lat_index, lon_index_west:lon_index] * np.sign(lat_field[lat_index])
    v_vel_east  = v_vel_field[lat_index, lon_index+1:lon_index_east+1] * np.sign(lat_field[lat_index])

    #Minimum acquired wind speed
    vel_crit        = 7.5 

    return ((np.min(u_vel_north) > -vel_crit) or (np.max(u_vel_south) < vel_crit) 
            or (np.max(v_vel_east) < vel_crit) or (np.min(v_vel_west) > -vel_crit))


def distance_adjacent_RV_maxima(PSL_min_lon, PSL_min_lat, lon, lat, track):
    """calculate distance to RV maxima from (extrapolated) last point of track"""
    
    if track.ndim > 1:  # track has at least two timestamps
        lon_1, lon_2    = track[-1, 1], track[-2, 1]
        lat_1, lat_2    = track[-1, 2], track[-2, 2]
        
        # linearly extrapolate to new point (x)
        lon_x   = lon_1 + (lon_1 - lon_2)  
        lat_x   = lat_1 + (lat_1 - lat_2)

        # determine distance from (x) to each new candidate point (p)
        distance    = np.zeros(len(PSL_min_lat)) + 1000.0
        for min_i in range(len(PSL_min_lat)):
            lon_p   = lon[PSL_min_lon[min_i]]
            lat_p   = lat[PSL_min_lat[min_i]]
            distance[min_i] = Distance(lon_x, lat_x, lon_p, lat_p) / 1000.0
            
    elif track.ndim == 1:  # track has only one timestamp
        lon_t   = track[1]
        lat_t   = track[2]

        # determine distance from last track point (t) to each new candidate point (p)
        distance    = np.zeros(len(PSL_min_lat)) + 1000.0
        for min_i in range(len(PSL_min_lat)):
            lon_p   = lon[PSL_min_lon[min_i]]
            lat_p   = lat[PSL_min_lat[min_i]]
            distance[min_i] = Distance(lon_t, lat_t, lon_p, lat_p) / 1000.0
    else:
        raise ValueError(f"cannot determine distance to adjacent maxima: \n{track.ndim=}\n{track=}")

    return distance


def fill_track_data(PSL_min_lon, PSL_min_lat, index_min, lon, lat, lat_weight, lon_ocn, lat_ocn, temp_anom, time_i, 
                    pres, pres_month, vor_850, U_10, SST_day, SST_month, u_vel_850, v_vel_850, u_vel_250, v_vel_250, weakening=False):
    #Get the indices and coordinates for the PSL minimum
    lon_index   = PSL_min_lon[index_min]
    lat_index   = PSL_min_lat[index_min]
    lon_pres    = lon[lon_index]
    lat_pres    = lat[lat_index]

    #Get the 850 hPa TEMP anom (8x8)
    temp_anom_TC    = temp_anom[index_min]

    #Get the pressure and vorticity
    pres_TC_day = pres[time_i, lat_index, lon_index]
    pres_TC_month   = pres_month[lat_index, lon_index]

    #Get the maximum vorticity near the PSL minima
    vor_850_TC_max  = PressureRVMaxima(lon, lat, lon_index, lat_index, vor_850[time_i])
    if weakening:
        vor_850_TC_max -= 1

    #Find the radius of maximum wind speeds
    vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i])

    if rad_vel_max > 400:
        #Reduce region of search
        vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i], reduced_radius = True)

    #Get the sea surface temperature
    SST_TC_day, SST_TC_month    = SeaSurfaceTemperatureTC(lon_ocn, lat_ocn, SST_day, SST_month, lon_pres, lat_pres)
    
    if SST_TC_day is np.ma.masked: 
        aboveland = True # TC candidate point above land
    else:
        aboveland = False

    #Determine vertical wind shear
    TC_shear    = VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i], u_vel_250[time_i], v_vel_250[time_i])

    #Find the maximum precipitation
    prec_max    = MaximumPrecipitation(lon, lat, lon_index, lat_index, prec[time_i])

    data_track  = np.ma.masked_all(13)
    data_track[0]   = time[time_i]  # time
    data_track[1]   = lon_pres      # longitude
    data_track[2]   = lat_pres      # latitude
    data_track[3]   = pres_TC_day   # pressure (3h)
    data_track[4]   = pres_TC_month # pressure (monthly)
    data_track[5]   = vel_max       # maximum 10m velocity
    data_track[6]   = rad_vel_max   # radius of maximum 10m velocity
    data_track[7]   = SST_TC_day    # sea surface temperature (daily)
    data_track[8]   = SST_TC_month  # sea surface temperature (monthly)
    data_track[9]   = vor_850_TC_max # maximum 850 hPa relative vorticity
    data_track[10]  = TC_shear      # vertical wind shear
    data_track[11]  = prec_max      # maximum precipitation
    data_track[12]  = temp_anom_TC  # 850 hPa core temperature anomaly
    
    return data_track, aboveland


def shorten_track_if_formed_above_land(track):
    """Remove parts of track such that the first 24 hours do not contain any land mask
    or only 24 hours of track remain (to be filtered out by track duration)
    """
    if len(track) <= 8: # in case of TC seeds
        if any(track[:,7].mask):
            return track[:0]
        else:
            return track

    for time_j in range(len(track) - 8):
        if any(track[time_j:time_j + 8, 7].mask):
            continue
        break  #No adjustments needed, first 24 hour do not contain land mask

    if time_j != 0:
        track   = track[time_j:]  #Adjust the track
        print(f'TRACK_ID_{track_i} = TRACK_ID_{track_i}[{time_j}:] (land mask detected)')
        exec(f'TRACK_ID_{track_i} = TRACK_ID_{track_i}[{time_j}:]', globals())
    
    return track


def has_24_hours_of_TC_strength(track):
    """returns if the track has 8 consecutive steps of TC-strength wind"""
    
    vel_max_counter = 0
                
    for time_j in range(len(track)):
        if track[time_j, 5] >= 17.0:  # at least 17 m/s 10m-wind
            vel_max_counter += 1
        else: 
            vel_max_counter = 0
        if vel_max_counter == 8:
            return True  # 8 consecutive time steps (24hr) of strong wind
    
    return False


#def shorten_track_until_TC_strength(track):
#    """end track before the first 24 consecutive hours of TC-strength wind (for TC seeds)"""
#
#    vel_max_counter = 0
#
#    print("shorten_track_until_TC_strength (threshold=17m/s)")
#    print("vel_max_counter, time_j, track[time_j,5]")
#    for time_j in range(len(track)):
#        if track[time_j, 5] >= 17.0:  # at least 17 m/s 10m-wind
#            vel_max_counter += 1
#        else:
#            vel_max_counter = 0
#        print(vel_max_counter, time_j, track[time_j,5])
#        if vel_max_counter == 8: # 8 consecutive time steps (24 hr) of strong wind
#            print(f"SHORTENED to: {track[:time_j-7,5]}")
#            return track[:time_j-7]
#    print("NO CHANGE")
#    return track


def shorten_track_until_TC_strength(track):
    """end track once V10 reaches 17 m/s"""
    if np.max(track[:,5]) > 17.0:
        idx = np.nonzero(track[:,5]>17.0)[0][0]
        track = track[:idx]
        print(f'TRACK_ID_{track_i} = TRACK_ID_{track_i}[:{idx}] (remove after v10>17m/s)')
        exec(f'TRACK_ID_{track_i} = TRACK_ID_{track_i}[:{idx}]', globals())
    return track


def has_24_hours_of_weak_shear(track):
    """returns whether the track has 8 consecutive steps with low windshear"""
    
    shear_counter   = 0

    for time_j in range(len(track)):
        if track[time_j, 10] <= 12.5:  # at most 12.5 m/s wind shear (250-850hPa)
            shear_counter += 1
        else:
            shear_counter = 0

        if shear_counter == 8:
            return True  # 8 consecutive time steps (24hr) of low shear
    
    return False


def validate_TC(track):
    """Assess if finished track is a valid TC. For this, the following checks are performed:

    1) duration should be at least 48 hours
    2) after removing potential formation over land, the remaining duration should be at least 48 hours
    3) genesis location should be within [30S-30N]
    4) genesis SST should be >= 25C
    5) 10m wind should be >= 17m/s for at least 24 consecutive hours
    6) wind shear should be < 12.5 m/s for at least 24 consecutive hours

    Returns:
        bool : True if a valid TC
        reason : str, short description why it is valid or not
    """

    # remove tracks lasting less than 48 hours directly
    duration_track  = 3.0 if (track.ndim == 1) else (len(track) * 3.0)
    if duration_track < 48.0:
        return False, 'tooshortregular'

    # shorten track such that the first 24hrs do not contain any land mask
    # and repeat the above
    track = shorten_track_if_formed_above_land(track)
    duration_track  = len(track) * 3.0 # assuming duration always >= 6 hrs
    if duration_track < 48.0:
        return False, 'tooshortland'

    # genesis location should be in the tropics
    if np.abs(track[0, 2]) > 30.0:
        return False, 'origin'

    # genesis SST should be at least 25 degC
    if track[0, 7] < 25.0:
        return False, 'SST'

    # candidate must have at least 24 consecutive hours of TC strength
    if not has_24_hours_of_TC_strength(track):
        return False, 'velmax'

    # candidate must have at least 24 consecutive hours of low shear
    if not has_24_hours_of_weak_shear(track):
        return False, 'shear'

    # all checks are passed, this is a valid TC
    return True, 'valid'
    
    
def validate_seed(track):
    """Assess if finished track is a valid TC seed. For this, the following checks are performed:

    1) duration should be at least 12 hours
    2) after removing potential formation over land, the remaining duration should be at least 12 hours
    3) after removing track data once V10 has reached 17 m/s, the remaining duration should be at least 12 hours
    4) genesis location should be within [30S-30N]
    5) genesis SST should be >= 25C

    Returns:
        bool : True if a valid TC
        reason : str, short description why it is valid or not
    """

    # remove tracks lasting less than 48 hours directly
    duration_track  = 3.0 if (track.ndim == 1) else (len(track) * 3.0)
    if duration_track < 12.0:
        return False, 'tooshortregular'

    # shorten track such that the first 24hrs do not contain any land mask
    # and repeat the above
    track = shorten_track_if_formed_above_land(track)
    duration_track  = len(track) * 3.0 # assuming duration always >= 6 hrs
    if duration_track < 12.0:
        return False, 'tooshortland'

    track = shorten_track_until_TC_strength(track)
    duration_track = len(track) * 3.0
    if duration_track < 12.0:
        return False, 'tooshortintensity'

    # genesis location should be in the tropics
    if np.abs(track[0, 2]) > 30.0:
        return False, 'origin'

    # genesis SST should be at least 25 degC
    if track[0, 7] < 25.0:
        return False, 'SST'

    ## candidate must have at least 24 consecutive hours of TC strength
    #if not has_24_hours_of_TC_strength(track):
    #    return False, 'velmax'

    ## candidate must have at least 24 consecutive hours of low shear
    #if not has_24_hours_of_weak_shear(track):
    #    return False, 'shear'

    # all checks are passed, this is a valid TC
    return True, 'valid'


def validate(track):
    """check if track is valid 
    (for more info, see subfunctions)
    """
    if SEEDS:
        isvalid, reason = validate_seed(track)
    else:
        isvalid, reason = validate_TC(track)
    
    return isvalid, reason


def remove_candidate(track_i, track_ID_active, reason):
    """remove candidate from active tracks"""
    track_ID_active.remove(track_i)
    print(f'del TRACK_ID_{track_i} ({reason})')
    exec(f'del TRACK_ID_{track_i}', globals()) # free up memory
    exec(f'counter_{reason} += 1', globals()) # diagnostic counter 
    return

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------  

print(f"{len(files)=}")
print(files[0])
print(files[-1])

# output file
print(f"output file: {outdir}/{outfile}")
if os.path.exists(f'{outdir}/{outfile}'):
    raise ValueError(f'{outfile} already exists, quitting...')

print(f"Track TC seeds: {SEEDS}")

tmpdir = os.path.join(os.environ['TMPDIR'], f'{experiment}.{run_ensemble}')
print(f"{tmpdir=}")
if not os.path.exists(tmpdir):
    os.mkdir(tmpdir)

#-----------------------------------------------------------------------------------------

with netcdf.Dataset(files[0],'r') as fh:
    time_units = fh['time'].units
    time_calendar = fh['time'].calendar
    time_slice = [cftime.datetime.strptime(t, "%Y-%m-%d", time_calendar) 
        for t in time_slice]

with netcdf.Dataset(gridfile, 'r') as fh:
    lon_grid = fh.variables['lon'][:]        # Longitude
    lat_grid = fh.variables['lat'][126:642]  # Latitude   
    grid_x   = fh.variables['DX'][126:642]   # Zonal length of grid cell  
    grid_y   = fh.variables['DY'][126:642]   # Meridional length of grid cell     

lon_2, grid_x   = PeriodicBoundaries2D(lon_grid, lat_grid, grid_x)
lon, grid_y = PeriodicBoundaries2D(lon_grid, lat_grid, grid_y)

#Generate for each grid cell the difference length
x_diff      = 0.5 * grid_x[:,2:] + 0.5 * grid_x[:,:-2] + grid_x[:,1:-1]
y_diff      = 0.5 * grid_y[2:,:] + 0.5 * grid_y[:-2,:] + grid_y[1:-1,:]

#Last elements can not be determined
x_diff  = x_diff[1:-1,:]
y_diff  = y_diff[:,1:-1]

#-----------------------------------------------------------------------------------------

track_ID    = 0           # counter for track ID
track_ID_active = []      # currently active track IDs
track_ID_list   = []      # valid track IDs (which will be saved)
track_attrs = ['']*13     # hold nc attributes for all 13 stored variables 

day_counter = 0           # keep track of total number of processed days
current_month   = 13      # Keep track of current month (for SSTs)
tasks = chunk_data(files) # daily chunks for all times in CAM data

# diagnostic counters
counter_notupdated = 0    # no new data added (could not match with new maximum or used optional maxima twice) 
counter_tooshortregular=0 # less than 48 hours
counter_tooshortland = 0  # less than 48 hours after removing part above land
counter_tooshortintensity=0 # less than 12 hours after removing TC-intensity part (only seeds)
counter_origin = 0        # no tropical origin
counter_velmax = 0        # no 24hr (consecutively) of at least 17m/s wind
counter_shear = 0         # no 24hr (consecutively) with less than 12.5m/s shear
counter_SST = 0           # SST less than 25C at genesis

#-----------------------------------------------------------------------------------------

# loop through all days in CAM
for task in tasks:        
    checkdate, filelist, t0, t1 = task
    
    # skip CAM time chunks outside time_slice
    checkdatetime = cftime.datetime.strptime(checkdate, "%Y-%m-%d", time_calendar)
    if checkdatetime < time_slice[0] or checkdatetime >= time_slice[-1]:
        print(f'skipping {checkdate}, out of time_slice')
        continue
        
    # match chunk with RV_Max file
    RVfile = f'{directory}/RV_Max/RV_Max_Coordinates_{checkdate}.nc'
    if not os.path.exists(RVfile):
        print(f"WARNING: skipping {checkdate}, no RVmax file")
        continue
    try:
        filenames = [files[fid] for fid in filelist]
    except IndexError:
        print(f"WARNING: task {task} failed. {len(files)=}")
        continue
        
    # start processing / read CAM data for a specific day
    day_counter += 1
    (time, timestamp, date, lon, lat, lat_weight, pres, U_10, vor_850, u_vel_850, v_vel_850, 
        u_vel_250, v_vel_250, prec) = ReadinData(filenames, x_diff, y_diff, t0, t1)

    # loop through CAM timestamps for this day
    for time_i in range(len(time)):
        print(f"\n{time_i} {timestamp[time_i]} ({time[time_i]} {time_units})")
        assert date[time_i] == checkdate, f"{date[time_i]=} not equal to {checkdate=}"
        
        #---------------------------------------------------------------------------------------
        #  read data
        #---------------------------------------------------------------------------------------
        
        # read RV_Max file
        try:
            PSL_min_lon, PSL_min_lat, temp_anom = ReadinData_RV_Max(RVfile, timestamp[time_i])
        except FileNotFoundError:
            print(f"WARNING: failed to read {RVfile}.")
            continue

        # read PSLmon, SSTmon and SST data if entering a new month
        time_year, time_month, time_day = (int(t) for t in date[time_i].split("-"))
        if time_month != current_month:
            pres_month = ReadinDataPSLMonth(time_month, time_year)
            time_ocn, lon_ocn, lat_ocn, SST, SST_month  = ReadinDataSST(time_month, time_year)
            current_month = time_month
        
        # select SSTday from SST
        if time_ocn is not None: 
            time_index  = np.argmin(np.abs(time_ocn-timestamp[time_i]))
            time_diff   = abs(time_ocn[time_index]-timestamp[time_i]).total_seconds()/3600
            if  time_diff > 13:
                print(f"WARNING: {time_diff} hr difference between current time {timestamp[time_i]}"
                    + f" and ocean time stamp {time_ocn[time_index]}.\n{time_ocn=}")
            SST_day     = SST[time_index]
        else:
            SST_day = None
            
        #---------------------------------------------------------------------------------------
        #  mark RV maxima beyond 60N/S or without closed circulation as 'optional' (*_opt)
        #---------------------------------------------------------------------------------------
        
        remove_index    = []  # indices of candidates that will be removed (in PSL_min_lon/lat)
        counter_bounds = 0    # only used for printing diagnostics
        counter_velocity = 0  # only used for printing diagnostics

        for min_i in range(len(PSL_min_lat)):   #Get the indices for the RV maxima
            lat_index   = PSL_min_lat[min_i]
            lon_index   = PSL_min_lon[min_i]

            if out_of_latbounds(lat_index, lat): # Low at 60S or 60N
                counter_bounds += 1
                remove_index.append(min_i)
                continue
    
            # storm-relative 10m tangential (cyclonic) wind speed lower than 7.5 m/s (0 for seeds)
            if has_no_closed_circulation(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i]):
                counter_velocity += 1
                remove_index.append(min_i)
                continue

        # Get the indices which are still optional
        optional_index  = [i for i in range(len(PSL_min_lon)) if i not in remove_index]
        PSL_min_lon_opt = np.delete(PSL_min_lon, optional_index)
        PSL_min_lat_opt = np.delete(PSL_min_lat, optional_index)
        temp_anom_opt   = np.delete(temp_anom, optional_index)

        # Remove the false RV maxima from array
        PSL_min_lon     = np.delete(PSL_min_lon, remove_index)
        PSL_min_lat     = np.delete(PSL_min_lat, remove_index)
        temp_anom   = np.delete(temp_anom, remove_index)
        
        print(f"RVfile[{time_i}]: {len(PSL_min_lon)} valid and {len(PSL_min_lon_opt)} optional ({counter_bounds} @60N/S and {counter_velocity} weak (<7.5m/s))")
        
        #-----------------------------------------------------------------------------------------
        # Loop through active tracks having at least two time steps and try to match with
        # regular or optional RV maxima.
        #-----------------------------------------------------------------------------------------
        
        #print(f"{len(track_ID_active)} active tracks from previous step, matching with new maxima...")
        counter_tooshort = 0
        counter_validmatch = 0
        counter_optmatch = 0
        counter_secondweakening = 0
        
        for track_i in track_ID_active: # loop over active tracks (from previous data chunk)
            exec(f'track = TRACK_ID_{track_i}')
            
            if track.ndim == 1:  # only one time stamp, continue (see next step)
                counter_tooshort += 1
                continue
                
            distance = distance_adjacent_RV_maxima(PSL_min_lon, PSL_min_lat, lon, lat, track)
            
            if np.any(distance <= 200) == True:  # Match with regular RV maximum
                counter_validmatch += 1

                #Check the previous RV at 850 (in case of temporarily weakening)
                vor_850_TC_prev = track[-1, 9]

                if vor_850_TC_prev < 0:
                    #Negative vorticity implies temporarily weakening, but now a valid PSL minima is found
                    #set previous RV to correct value (raise by 1) and continue
                    print(f'unflag TRACK_ID_{track_i} (remove weakening flag)')
                    exec('TRACK_ID_'+str(track_i)+'[-1, 9] += 1')

                #Check for the closest PSL minima in the neighbourhood of extrapolated PSL minima
                index_min   = np.argmin(distance)
                
                data_track,_ = fill_track_data(PSL_min_lon, PSL_min_lat, index_min, lon, lat, lat_weight, lon_ocn, lat_ocn, temp_anom, time_i, 
                    pres, pres_month, vor_850, U_10, SST_day, SST_month, u_vel_850, v_vel_850, u_vel_250, v_vel_250)

                #Get the correct ID
                #print(f'app TRACK_ID_{track_i} = {np.array_str(data_track[1:3],precision=2)}')
                exec(f'TRACK_ID_{track_i} = np.ma.vstack([TRACK_ID_{track_i}, data_track])')

                #Remove the already assigned PSL minima from array
                PSL_min_lon     = np.delete(PSL_min_lon, index_min)
                PSL_min_lat     = np.delete(PSL_min_lat, index_min)
                temp_anom   = np.delete(temp_anom, index_min)

            else:  # No match with regular RV maximum, check if previous RV maximum was optional
                
                vor_850_TC_prev = track[-1, 9]
                if vor_850_TC_prev < 0:  # Negative RV850 implies the previous point was optional -> end of track
                    counter_secondweakening += 1
                    # do no update TC track and remove previous point
                    print(f'end TRACK_ID_{track_i} (second weakening)')
                    exec(f'TRACK_ID_{track_i} = TRACK_ID_{track_i}[:-1]')

                else:  # Previous value was a regular part of the track, try to match with optional RV maxima
                    distance = distance_adjacent_RV_maxima(PSL_min_lon_opt, PSL_min_lat_opt, lon, lat, track)
                    
                    if np.any(distance <= 200) == True: # matched with new (optional) RV maximum
                        counter_optmatch += 1
                        index_min   = np.argmin(distance)
                        data_track,_ = fill_track_data(PSL_min_lon_opt, PSL_min_lat_opt, index_min, lon, lat, lat_weight, lon_ocn, lat_ocn, temp_anom_opt, time_i, 
                            pres, pres_month, vor_850, U_10, SST_day, SST_month, u_vel_850, v_vel_850, u_vel_250, v_vel_250, weakening=True)

                        # append new data
                        print(f'flag TRACK_ID_{track_i} = {np.array_str(data_track[1:3],precision=2)} (first weakening)')
                        exec(f'TRACK_ID_{track_i} = np.ma.vstack([TRACK_ID_{track_i}, data_track])')

                        #Remove the already assigned PSL minima from array
                        PSL_min_lon_opt = np.delete(PSL_min_lon_opt, index_min)
                        PSL_min_lat_opt = np.delete(PSL_min_lat_opt, index_min)
                        temp_anom_opt   = np.delete(temp_anom_opt, index_min)

        
        #print(f"{len(track_ID_active)-counter_tooshort} active tracks (>2) matched {counter_validmatch} valid, {counter_optmatch} optional ({counter_secondweakening} weak removed)")
        
        #-----------------------------------------------------------------------------------------
        # Loop through active tracks having only one time steps and try to match with
        # regular RV maxima.
        #-----------------------------------------------------------------------------------------
        
        remove_index    = []

        for track_i in track_ID_active:
            exec(f'track = TRACK_ID_{track_i}')

            if track.ndim > 1:
                continue
            
            distance = distance_adjacent_RV_maxima(PSL_min_lon, PSL_min_lat, lon, lat, track)
            
            if np.any(distance <= 200) == True:  # found regular RV maximum that matches track
                index_min   = np.argmin(distance)
                data_track,_ = fill_track_data(PSL_min_lon, PSL_min_lat, index_min, lon, lat, lat_weight, lon_ocn, lat_ocn, temp_anom, time_i, 
                    pres, pres_month, vor_850, U_10, SST_day, SST_month, u_vel_850, v_vel_850, u_vel_250, v_vel_250)

                #print(f'app TRACK_ID_{track_i} = {np.array_str(data_track[1:3],precision=2)}')
                exec(f'TRACK_ID_{track_i} = np.ma.vstack([TRACK_ID_{track_i}, data_track])')
                remove_index.append(index_min)  #Current RV maxima will not start a new track
                
        PSL_min_lon     = np.delete(PSL_min_lon, remove_index)
        PSL_min_lat     = np.delete(PSL_min_lat, remove_index)
        temp_anom   = np.delete(temp_anom, remove_index)
        
        #print(f"{counter_tooshort} active tracks (1) matched {len(remove_index)} valid")
        
        #-----------------------------------------------------------------------------------------
        #  The remaining regular RV maxima will start new tracks
        #-----------------------------------------------------------------------------------------
        
        for min_i in range(len(PSL_min_lat)):
            
            data_track, aboveland = fill_track_data(PSL_min_lon, PSL_min_lat, min_i, lon, lat, lat_weight, lon_ocn, lat_ocn, temp_anom, time_i, 
                    pres, pres_month, vor_850, U_10, SST_day, SST_month, u_vel_850, v_vel_850, u_vel_250, v_vel_250)
            
            if aboveland:
                continue  # TC forms above land, discard

            #Initiate new track for new low
            print(f'ini TRACK_ID_{track_ID} = {np.array_str(data_track[1:3],precision=2)}')
            exec(f'TRACK_ID_{track_ID} = data_track')

            #Save all the ID's, some ID's will be removed (too short track etc.)
            track_ID_active.append(track_ID)
            
            #Raise track ID
            track_ID    += 1

        #print(f"{len(PSL_min_lon)} new active tracks (1) created -> {len(track_ID_active)=}.")
        
        #-----------------------------------------------------------------------------------------
        # Now going through all recently finished tracks (i.e. no new data added), perform some 
        # final checks and either archive or delete them.
        #-----------------------------------------------------------------------------------------
        
        # diagnostic counters
        #counter_notupdated = 0
        #counter_tooshortregular = 0
        #counter_tooshortland = 0
        #counter_origin = 0
        #counter_velmax = 0
        #counter_shear = 0
        #counter_SST = 0
    
        for track_i in track_ID_active:
            exec(f'track = TRACK_ID_{track_i}')

            # check if track is still active
            time_end_track = track[0] if (track.ndim == 1) else track[-1,0]
            if time[time_i] == time_end_track:
                continue  # Still an active track, keep active
            
            counter_notupdated += 1
             
            isvalid, msg = validate(track)
            if isvalid:
                print(f"save TRACK_ID_{track_i} ({msg})")
                track_ID_list.append(track_i)
                track_ID_active.remove(track_i)
            else:
                remove_candidate(track_i, track_ID_active, msg)
            
            ## remove tracks lasting less than 48 hours directly
            #duration_track  = 3.0 if (track.ndim == 1) else (len(track) * 3.0)
            #if duration_track < 48.0:
            #    #counter_tooshortregular += 1
            #    remove_candidate(track_i, track_ID_active, 'tooshortregular')
            #    continue

            #
            ##print(f"{track_i=} finished: {duration_track=}hr")
            ##Storm lasts more than 2 days

            #track = shorten_track_if_formed_above_land(track)
            #
            #duration_track  = len(track) * 3.0
            #if duration_track < 48.0:
            #    #counter_tooshortland += 1
            #    remove_candidate(track_i, track_ID_active, 'tooshortland')
            #    continue

            ## genesis location should be in the tropics 
            #if np.abs(track[0, 2]) > 30.0:
            #    #counter_origin += 1
            #    remove_candidate(track_i, track_ID_active, 'origin')
            #    continue

            ## genesis SST should be at least 25 degC
            #if track[0, 7] < 25.0:
            #    #counter_SST += 1
            #    remove_candidate(track_i, track_ID_active, 'SST')
            #    continue

            ## candidate must have at least 24 consecutive hours of TC strength
            #if not has_24_hours_of_TC_strength(track):
            #    #counter_velmax += 1
            #    remove_candidate(track_i, track_ID_active, 'velmax')
            #    continue

            ## candidate must have at least 24 consecutive hours of low shear
            #if not has_24_hours_of_weak_shear(track):
            #    #counter_shear += 1
            #    remove_candidate(track_i, track_ID_active, 'shear')
            #    continue

            ## if all checks are passed, add to archive list
            #print(f"store TRACK_ID_{track_i} (valid track, {duration_track=}hr, archiving...)")
            #track_ID_list.append(track_i)
            #track_ID_active.remove(track_i)

        print(f"{len(track_ID_active)} active tracks")
        print(f"{len(track_ID_list)} archived tracks (valid)")



    print(f"\n{counter_notupdated} tracks finished -> {len(track_ID_active)=}\n  discard: {counter_tooshortregular+counter_tooshortland+counter_tooshortintensity} short "
        + f"({counter_tooshortregular} reg. {counter_tooshortland} land, {counter_tooshortintensity} intensity), {counter_origin} non trop., "
        + f"{counter_velmax} weak, {counter_shear} shear, {counter_SST} SST")

#-----------------------------------------------------------------------------------------
#  Now repeat for tracks that are still active at the end of the run
#-----------------------------------------------------------------------------------------

#counter_tooshortregular = 0
#counter_tooshortland = 0
#counter_origin = 0
#counter_velmax = 0
#counter_shear = 0
#counter_SST = 0  

for track_i in track_ID_active:
    exec(f'track = TRACK_ID_{track_i}')

    isvalid, msg = validate(track)
    if isvalid:
        print(f"save TRACK_ID_{track_i} ({msg})")
        track_ID_list.append(track_i)
    else:
        remove_candidate(track_i, track_ID_active, msg)

    ## remove tracks lasting less than 48 hours directly
    #duration_track  = 3.0 if (track.ndim == 1) else (len(track) * 3.0)
    #if duration_track < 48.0:
    #    #counter_tooshortregular += 1
    #    remove_candidate(track_i, track_ID_active, 'tooshortregular')
    #    continue

    ##print(f"{track_i=} finished: {duration_track=}hr")
    ##Storm lasts more than 2 days

    #track = shorten_track_if_formed_above_land(track)

    #duration_track  = len(track) * 3.0
    #if duration_track < 48.0:
    #    #counter_tooshortland += 1
    #    remove_candidate(track_i, track_ID_active, 'tooshortland')
    #    continue

    ## genesis location should be in the tropics 
    #if np.abs(track[0, 2]) > 30.0:
    #    #counter_origin += 1
    #    remove_candidate(track_i, track_ID_active, 'origin')
    #    continue

    ## genesis SST should be at least 25 degC
    #if track[0, 7] < 25.0:
    #    #counter_SST += 1
    #    remove_candidate(track_i, track_ID_active, 'SST')
    #    continue

    ## candidate must have at least 24 consecutive hours of TC strength
    #if not has_24_hours_of_TC_strength(track):
    #    #counter_velmax += 1
    #    remove_candidate(track_i, track_ID_active, 'velmax')
    #    continue

    ## candidate must have at least 24 consecutive hours of low shear
    #if not has_24_hours_of_weak_shear(track):
    #    #counter_shear += 1
    #    remove_candidate(track_i, track_ID_active, 'shear')
    #    continue

    ## if all checks are passed, add to archive list
    #print(f"save TRACK_ID_{track_i} (valid), {duration_track=}hr, archiving...")
    #track_ID_list.append(track_i)
    
print(f"Tracks still active at the end ({len(track_ID_active)}):\n  {counter_tooshortregular+counter_tooshortland+counter_tooshortintensity} short "
    + f"({counter_tooshortregular} reg. {counter_tooshortland} land, {counter_tooshortintensity} intensity), {counter_origin} non trop., "
    + f"{counter_velmax} weak, {counter_shear} shear, {counter_SST} SST")       
print(f"Total number of valid stored tracks: {len(track_ID_list)}")



#-----------------------------------------------------------------------------------------
#  Merge all valid tracks
#-----------------------------------------------------------------------------------------

#Determine the maximum track time of all TCs
time_max    = 0

for track_i in track_ID_list:
    exec(f'track = TRACK_ID_{track_i}')
    if len(track) > time_max:
        time_max    = len(track)
        
#Raise time max by 1 (easy to remove masked objects for post-processing)
time_max    += 1

# all track data
track_all   = np.ma.masked_all((len(track_ID_list), time_max, 13))

for track_i in range(len(track_ID_list)):
    #Save data to general array
    exec('track = TRACK_ID_'+str(track_ID_list[track_i]))
    print(f'saving TRACK_ID_{track_i}: {len(track)=}')
    track_all[track_i, :len(track)] = track
    
track_all[:,:,7:9][track_all[:,:,7:9]==100] = np.nan # missing SST/SSTmon values

#-----------------------------------------------------------------------------------------
#  Save to NetCDF
#-----------------------------------------------------------------------------------------

HEAT_data = netcdf.Dataset(f'{outdir}/{outfile}', 'w')

HEAT_data.createDimension('id', len(track_ID_list))
HEAT_data.createDimension('dtime', time_max)
HEAT_data.createDimension('data', 13)

HEAT_data.createVariable('id', float, ('id'), zlib=True)
HEAT_data.createVariable('dtime', float, ('dtime'), zlib=True)
HEAT_data.createVariable('data', float, ('data'), zlib=True)
HEAT_data.createVariable('TC_tracks', float, ('id', 'dtime', 'data'), 
                         fill_value=9.969209968386869e+36, zlib=True)
HEAT_data.createVariable('num_days', int, zlib=True)

HEAT_data.variables['id'].long_name = 'Track number'
HEAT_data.variables['dtime'].long_name  = 'Lifetime TC'
HEAT_data.variables['dtime'].units      = 'days'
HEAT_data.variables['data'].description = track_attrs
HEAT_data.variables['num_days'].long_name = 'Total number of analysed days'

#Writing data to correct variable   
HEAT_data.variables['id'][:]    = np.arange(len(track_ID_list)) + 1
HEAT_data.variables['dtime'][:]     = np.arange(time_max) * 0.125
HEAT_data.variables['data'][:]      = np.arange(13)
HEAT_data.variables['TC_tracks'][:] = track_all
HEAT_data.variables['num_days'][:]  = day_counter 

HEAT_data.close()
print(f"Results saved to {outdir}/{outfile}")

#-----------------------------------------------------------------------------------------
#  Delete temporary files
#-----------------------------------------------------------------------------------------
for file in sorted(glob.glob(os.path.join(tmpdir,'*'))):
    os.remove(file)

os.rmdir(tmpdir)
