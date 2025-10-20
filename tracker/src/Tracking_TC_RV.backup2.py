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
    RVSEARCHRADIUS = 200 # (km) search radius for TC at previous time step
    NOUT            = 8 # number of time steps per day
    time_slice      = ('2003-01-01','2008-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.nc'
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
    NOUT            = 8 # number of time steps per day
    time_slice      = ('2003-01-01','2008-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.nc'
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
    RVSEARCHRADIUS = 200 # (km) search radius for TC at previous time step
    NOUT            = 8 # number of time steps per day
    time_slice      = ('2093-01-01','2098-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.nc'
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
    NOUT            = 8 # number of time steps per day
    time_slice      = ('2093-01-01','2098-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.nc'
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
    NOUT            = 8 # number of time steps per day
    time_slice      = ('2093-01-01','2098-01-01') # time period to analyse
    outfile         = f'TC_tracker_results.{experiment}.started_{run_year_start}.{run_ensemble:03d}.nc'
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
    return chunks
    
    
def getattstr(dataset, name):
    ncvar = dataset[name]
    attrs = {a:getattr(ncvar,a) for a in ncvar.ncattrs()}
    if name == 'PSL':
        attrs['units'] = 'h'+attrs['units']
        if 'cell_methods' in attrs and 'time: mean' in attrs['cell_methods']:
            name = 'PSLmon' 
    return f"{name}:{attrs}"


def ReadinData(filenames, x_diff, y_diff, t1, t2):
    dt = 1/NOUT # output time step in days
    with netcdf.MFDataset(filenames, 'r', aggdim='time') as fh:
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
    
    
def ReadinDataPRECT(timestamp):
    """Get precipitation from separate file stream, linearly interpolating to timestamp"""
    files = sorted(glob.glob(f'{datadir["atm"]}/{experiment_name}.cam2.{stream_prec}.*.nc'))
    dates = [file[-19:-3] for file in files]
    date0, date1 = timestamp[0].strftime('%Y-%m-%d'), timestamp[-1].strftime('%Y-%m-%d')
    fid0 = max(np.searchsorted(dates, date0)-1,0)
    fid1 = np.searchsorted(dates, date1)+1
    with netcdf.MFDataset(files[fid0:fid1], 'r') as fh: 
        nc_time = fh['time'][:] # Time (end of interval)
        nc_tbnd = fh['time_bnds'][:] # time bounds
        nc_prec = fh['PRECT'][:,127:641] # precipitation
        if 'time: mean' in getattr(fh['PRECT'], 'cell_methods', ''):
            nc_time = nc_tbnd.mean(axis=-1) # Center time if data represents time average
            if np.diff(nc_tbnd[0]) == 0.125: # average 3hrly data to 6hrly
                nc_time = np.mean([nc_time[1::2],nc_time[:-1:2]], axis=0)
                nc_prec = np.mean([nc_prec[1::2],nc_prec[:-1:2]], axis=0)
        times = cftime.date2num(timestamp, fh['time'].units, fh['time'].calendar) # convert timestamp to same units
        prec = interp1d(nc_time, nc_prec, axis=0, bounds_error=False)(times)
        if times[0] < nc_time[0] or times[-1] > nc_time[-1]:
            ids_invalid = (times < nc_time[0]) | (times > nc_time[-1])
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

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------  

# input files
files = sorted(glob.glob(f'{datadir["atm"]}/{experiment_name}.cam2.{stream}.*.nc'))
# files = files[6:373]
#files = files[:7]; print(f'line 598 (testrun): using {len(files)} files')
#files = files[-2:]; print(f'line 599 (testrun): using {len(files)} files')
#files = files[-5:]; print(f'line 600 (testrun): using {len(files)} files')

print(f"{len(files)=}")
print(files[0])
print(files[-1])

# output file
print(f"output file: {directory}/{outfile}")
if os.path.exists(f'{directory}/{outfile}'):
    raise ValueError(f'{outfile} already exists, quitting...')

#-----------------------------------------------------------------------------------------

with netcdf.Dataset(files[0],'r') as fh:
    time_units = fh['time'].units
    time_calendar = fh['time'].calendar
    time_slice = [cftime.datetime.strptime(t, "%Y-%m-%d", time_calendar) 
        for t in time_slice]

with netcdf.Dataset(gridfile, 'r') as fh:
    #Writing data to correct variable
    lon_grid    = fh.variables['lon'][:]    #Longitude
    lat_grid    = fh.variables['lat'][126:642]  #Latitude   
    grid_x      = fh.variables['DX'][126:642]   #Zonal length of grid cell  
    grid_y      = fh.variables['DY'][126:642]   #Meridional length of grid cell     

lon_2, grid_x   = PeriodicBoundaries2D(lon_grid, lat_grid, grid_x)
lon, grid_y = PeriodicBoundaries2D(lon_grid, lat_grid, grid_y)

#Generate for each grid cell the difference length
x_diff      = 0.5 * grid_x[:, 2:] + 0.5 * grid_x[:, :-2] + grid_x[:, 1:-1]
y_diff      = 0.5 * grid_y[2:] + 0.5 * grid_y[:-2] + grid_y[1:-1]

#Last elements can not be determined
x_diff  = x_diff[1:-1]
y_diff  = y_diff[:, 1:-1]

#Counter for track ID
track_ID    = 0
track_ID_active = []    #Currently active track IDs
track_ID_list   = []    #Valid track IDs

track_attrs = ['']*13 #np.empty(13, dtype=str) # 

#-----------------------------------------------------------------------------------------

day_counter = 0 # keep track of total number of processed days
current_month   = 13 # Keep track of current month (for SSTs)
tasks = chunk_data(files) # daily chunks
for task in tasks:
    checkdate, filelist, t0, t1 = task

    checkdatetime = cftime.datetime.strptime(checkdate, "%Y-%m-%d", time_calendar)
    if checkdatetime < time_slice[0] or checkdatetime >= time_slice[-1]:
        print(f'skipping {checkdate}, out of time_slice')
        continue
    RVfile = f'{directory}/RV_Max/RV_Max_Coordinates_{checkdate}.nc'
    if not os.path.exists(RVfile):
        print(f"WARNING: skipping {checkdate}, no RVmax file")
        continue
    try:
        filenames = [files[fid] for fid in filelist]
    except IndexError:
        print(f"WARNING: task {task} failed. {len(files)=}")
        continue
    day_counter += 1
        
    #Get the fields
    #print(f"{[os.path.basename(f) for f in filenames]} steps {t0}-{t1}") # this line CAUSES ERROR!
    (time, timestamp, date, lon, lat, lat_weight, pres, U_10, vor_850, u_vel_850, v_vel_850, 
        u_vel_250, v_vel_250, prec) = ReadinData(filenames, x_diff, y_diff, t0, t1)

    for time_i in range(len(time)):
        #loop over each 3-hourly field
        print(f"\n{time_i} {timestamp[time_i]}")
        assert date[time_i] == checkdate, f"{date[time_i]=} not equal to {checkdate=}"
        
        #Get the coordinates for each PSL minima
        try:
            PSL_min_lon, PSL_min_lat, temp_anom = ReadinData_RV_Max(RVfile, timestamp[time_i])
        except FileNotFoundError:
            continue

        #Retain the daily and monthly-averaged SST
        time_year, time_month, time_day = (int(t) for t in date[time_i].split("-"))

        if time_month != current_month:
            #New month, get the corresponding sea-level pressure and SSTs
            pres_month                  = ReadinDataPSLMonth(time_month, time_year)
            time_ocn, lon_ocn, lat_ocn, SST, SST_month  = ReadinDataSST(time_month, time_year)
            current_month                   = time_month

        #Get the corresponding daily SST
        if time_ocn is not None:
            time_index  = np.argmin(np.abs(time_ocn-timestamp[time_i]))
            time_diff   = abs(time_ocn[time_index]-timestamp[time_i]).total_seconds()/3600
            if  time_diff > 13:
                print(f"WARNING: {time_diff} hr difference between current time {timestamp[time_i]}"
                    + f" and ocean time stamp {time_ocn[time_index]}.\n{time_ocn=}")
            SST_day     = SST[time_index]
        else:
            SST_day = None
        
        #Keep track of RV maxima with a sea-level pressure minima in the neighbourhood
        remove_index    = []
        counter_bounds = 0
        counter_velocity = 0
        for min_i in range(len(PSL_min_lat)):   #Get the indices for the RV maxima
            lat_index   = PSL_min_lat[min_i]
            lon_index   = PSL_min_lon[min_i]

            if lat_index == 0 or lat_index == len(lat) - 1:
                #Low at 60S or 60N
                counter_bounds += 1
                remove_index.append(min_i)
                continue
    
            #Get a part of the field
            lon_field, lat_field, weight, u_vel_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, u_vel_850[time_i])
            lon_field, lat_field, weight, v_vel_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, v_vel_850[time_i])

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

            #Retain the velocity (add -1 for SH)
            u_vel_north = u_vel_field[lat_index+1:lat_index_north+1, lon_index] * np.sign(lat_field[lat_index])
            u_vel_south = u_vel_field[lat_index_south:lat_index, lon_index] * np.sign(lat_field[lat_index])
            v_vel_west  = v_vel_field[lat_index, lon_index_west:lon_index] * np.sign(lat_field[lat_index])
            v_vel_east  = v_vel_field[lat_index, lon_index+1:lon_index_east+1] * np.sign(lat_field[lat_index])

            #Minimum acquired wind speed
            vel_crit        = 7.5
    
            if np.min(u_vel_north) > -vel_crit or np.max(u_vel_south) < vel_crit or np.max(v_vel_east) < vel_crit or np.min(v_vel_west) > -vel_crit:
                #No closed cyclone circulation, discard RV maxima
                counter_velocity += 1
                remove_index.append(min_i)
                continue

        #-----------------------------------------------------------------------------------------

        #Get the indices which are still optional
        optional_index  = [i for i in range(len(PSL_min_lon)) if i not in remove_index]
        PSL_min_lon_opt = np.delete(PSL_min_lon, optional_index)
        PSL_min_lat_opt = np.delete(PSL_min_lat, optional_index)
        temp_anom_opt   = np.delete(temp_anom, optional_index)

        #Remove the false RV maxima from array
        PSL_min_lon     = np.delete(PSL_min_lon, remove_index)
        PSL_min_lat     = np.delete(PSL_min_lat, remove_index)
        temp_anom   = np.delete(temp_anom, remove_index)
        
        print(f"RVfile[{time_i}]: {len(PSL_min_lon)} valid and {len(PSL_min_lon_opt)} optional ({counter_bounds} @60N/S and {counter_velocity} weak (<7.5m/s))")
        
        #-----------------------------------------------------------------------------------------
        # Assign RV maxima to active tracks having at least two time steps.
        # On each time step, loop over active tracks (from previous chunk of data) and get last two 
        # coordinates. Extrapolate location and search for a new minimum within 200km.
        # If found, store TC data and remove that minimum from list
        # Else search in optional minimum list and
        print(f"{len(track_ID_active)} active tracks from previous step, matching with new maxima...")
        counter_tooshort = 0
        counter_validmatch = 0
        counter_optmatch = 0
        counter_secondweakening = 0
        
        for track_i in track_ID_active: # loop over active tracks (from previous data chunk)
            exec('track = TRACK_ID_'+str(track_i))
            
            if track.ndim == 1:
                #Only one time stamp, continue (see next step)
                counter_tooshort += 1
                continue

            #At least two time stamps, retain coordinates
            lon_1, lon_2    = track[-1, 1], track[-2, 1]
            lat_1, lat_2    = track[-1, 2], track[-2, 2]

            #Linearly extrapolate to new point (below for further processing)
            lon_1   = lon_1 + (lon_1 - lon_2)
            lat_1   = lat_1 + (lat_1 - lat_2)

            #Empty array for determining the distance
            distance    = np.zeros(len(PSL_min_lat)) + 1000.0

            for min_i in range(len(PSL_min_lat)):
    
                #Get the coordinates for potential RV maxima
                lon_2   = lon[PSL_min_lon[min_i]]
                lat_2   = lat[PSL_min_lat[min_i]]

                #Determine the distance in km
                distance[min_i] = Distance(lon_1, lat_1, lon_2, lat_2) / 1000.0

            if np.any(distance <= 200) == True:
                #There is a valid PSL minima in the neighbourhood
                counter_validmatch += 1

                #Check the previous RV at 850 (in case of temporarily weakening)
                vor_850_TC_prev = track[-1, 9]

                if vor_850_TC_prev < 0:
                    #Negative vorticity implies temporarily weakening, but now a valid PSL minima is found
                    #set previous RV to correct value (raise by 1) and continue
                    exec('TRACK_ID_'+str(track_i)+'[-1, 9] += 1')

                #Check for the closest PSL minima in the neighbourhood of extrapolated PSL minima
                index_min   = np.argmin(distance)

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
            
                #Find the radius of maximum wind speeds
                vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i])

                if rad_vel_max > 400:
                    #Reduce region of search
                    vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i], reduced_radius = True)

                #Get the sea surface temperature
                SST_TC_day, SST_TC_month    = SeaSurfaceTemperatureTC(lon_ocn, lat_ocn, SST_day, SST_month, lon_pres, lat_pres)

                #Determine vertical wind shear
                TC_shear    = VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i], u_vel_250[time_i], v_vel_250[time_i])

                #Find the maximum precipitation
                prec_max    = MaximumPrecipitation(lon, lat, lon_index, lat_index, prec[time_i])

                #Data for track (time, lon (pres), lat (pres), pres (3 h), pres (month), vel_max, rad. vel_max, SST (day), SST (month), vor max (1/s), TC shear (m/s), prec max (mm / day), 850 hpa temp anom (C))
                data_track  = np.ma.masked_all(13)
                data_track[0]   = time[time_i]
                data_track[1]   = lon_pres
                data_track[2]   = lat_pres
                data_track[3]   = pres_TC_day
                data_track[4]   = pres_TC_month
                data_track[5]   = vel_max
                data_track[6]   = rad_vel_max
                data_track[7]   = SST_TC_day
                data_track[8]   = SST_TC_month
                data_track[9]   = vor_850_TC_max
                data_track[10]  = TC_shear
                data_track[11]  = prec_max
                data_track[12]  = temp_anom_TC

                #Get the correct ID
                exec('TRACK_ID_'+str(track_i)+' = np.ma.vstack([TRACK_ID_'+str(track_i)+', data_track])')

                #Remove the already assigned PSL minima from array
                PSL_min_lon     = np.delete(PSL_min_lon, index_min)
                PSL_min_lat     = np.delete(PSL_min_lat, index_min)
                temp_anom   = np.delete(temp_anom, index_min)

            else:
                #No new valid point is found, use the optional PSL minima (only once)   
                #Check the previous RV at 850
                vor_850_TC_prev = track[-1, 9]

                if vor_850_TC_prev < 0:
                    counter_secondweakening += 1
                    #Negative vorticity implies that this is the second time in this loop, do no update TC track
                    #Remove the previous values as well
                    exec('TRACK_ID_'+str(track_i)+' = TRACK_ID_'+str(track_i)+'[:-1]')

                else:
                    #Previous value was still part of the track
                    #Search for local PSL minima from the optional ones
                    distance    = np.zeros(len(PSL_min_lat_opt)) + 1000.0

                    for min_i in range(len(PSL_min_lat_opt)):
                        #Get the coordinates for potential PSL minima
                        lon_2   = lon[PSL_min_lon_opt[min_i]]
                        lat_2   = lat[PSL_min_lat_opt[min_i]]

                        #Determine the distance in km
                        distance[min_i] = Distance(lon_1, lat_1, lon_2, lat_2) / 1000.0

                    if np.any(distance <= 200) == True:
                        counter_optmatch += 1
                        #Check for the closest PSL minima in the neighbourhood of extrapolated PSL minima
                        index_min   = np.argmin(distance)

                        #Get the indices and coordinates for the PSL minimum
                        lon_index   = PSL_min_lon_opt[index_min]
                        lat_index   = PSL_min_lat_opt[index_min]
                        lon_pres    = lon[lon_index]
                        lat_pres    = lat[lat_index]

                        #Get the 850 hPa TEMP anom (8x8)
                        temp_anom_TC    = temp_anom_opt[index_min]

                        #Get the pressure and vorticity
                        pres_TC_day = pres[time_i, lat_index, lon_index]
                        pres_TC_month   = pres_month[lat_index, lon_index]

                        #Get the maximum vorticity near the PSL minima
                        vor_850_TC_max  = PressureRVMaxima(lon, lat, lon_index, lat_index, vor_850[time_i])
                    
                        #Find the radius of maximum wind speeds
                        vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i])

                        if rad_vel_max > 400:
                            #Reduce region of search
                            vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i], reduced_radius = True)

                        #Get the sea surface temperature
                        SST_TC_day, SST_TC_month    = SeaSurfaceTemperatureTC(lon_ocn, lat_ocn, SST_day, SST_month, lon_pres, lat_pres)

                        #Determine vertical wind shear
                        TC_shear    = VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i], u_vel_250[time_i], v_vel_250[time_i])

                        #Find the maximum precipitation
                        prec_max    = MaximumPrecipitation(lon, lat, lon_index, lat_index, prec[time_i])

                        #Data for track (time, lon (pres), lat (pres), pres (3 h), pres (month), vel_max, rad. vel_max, SST (day), SST (month), vor max (1/s), TC shear (m/s), prec max (mm / day), 850 hpa temp anom (C))
                        data_track  = np.ma.masked_all(13)
                        data_track[0]   = time[time_i]
                        data_track[1]   = lon_pres
                        data_track[2]   = lat_pres
                        data_track[3]   = pres_TC_day
                        data_track[4]   = pres_TC_month
                        data_track[5]   = vel_max
                        data_track[6]   = rad_vel_max
                        data_track[7]   = SST_TC_day
                        data_track[8]   = SST_TC_month
                        data_track[9]   = vor_850_TC_max - 1.0  #Mark the vorticity as negative (temporarily)
                        data_track[10]  = TC_shear
                        data_track[11]  = prec_max
                        data_track[12]  = temp_anom_TC

                        #Get the correct ID
                        exec('TRACK_ID_'+str(track_i)+' = np.ma.vstack([TRACK_ID_'+str(track_i)+', data_track])')

                        #Remove the already assigned PSL minima from array
                        PSL_min_lon_opt = np.delete(PSL_min_lon_opt, index_min)
                        PSL_min_lat_opt = np.delete(PSL_min_lat_opt, index_min)
                        temp_anom_opt   = np.delete(temp_anom_opt, index_min)

        
        print(f"{len(track_ID_active)-counter_tooshort} active tracks (>2) matched {counter_validmatch} valid, {counter_optmatch} optional ({counter_secondweakening} weak removed)")
        
        #-----------------------------------------------------------------------------------------
        # Assign RV maxima to active tracks of length 1 
        # These will only be matched with valid indices.
        
        remove_index    = []

        for track_i in track_ID_active: # 
            #Loop over each active track
            exec('track = TRACK_ID_'+str(track_i))

            if track.ndim > 1:
                #More time shape, continue
                continue

            #Get the coordinates
            lon_1   = track[1]
            lat_1   = track[2]

            #Empty array for determining the distance
            distance    = np.zeros(len(PSL_min_lat)) + 1000.0
            
            for min_i in range(len(PSL_min_lat)):
                #Get the coordinates for the other potential points
                lon_2   = lon[PSL_min_lon[min_i]]
                lat_2   = lat[PSL_min_lat[min_i]]

                #Determine the distance in km
                distance[min_i] = Distance(lon_1, lat_1, lon_2, lat_2) / 1000.0

            if np.any(distance <= 200) == True:
                #Check for PSL minima in the neighbourhood
                index_min   = np.argmin(distance)

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

                #Find the radius of maximum wind speeds
                vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i])

                if rad_vel_max > 400:
                    #Reduce region of search
                    vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i], reduced_radius = True)

                #Get the sea surface temperature
                SST_TC_day, SST_TC_month    = SeaSurfaceTemperatureTC(lon_ocn, lat_ocn, SST_day, SST_month, lon_pres, lat_pres)

                #Determine vertical wind shear
                TC_shear    = VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i], u_vel_250[time_i], v_vel_250[time_i])
            
                #Find the maximum precipitation
                prec_max    = MaximumPrecipitation(lon, lat, lon_index, lat_index, prec[time_i])

                #Data for track (time, lon (pres), lat (pres), pres (3 h), pres (month), vel_max, rad. vel_max, SST (day), SST (month), vor max (1/s), TC shear (m/s), prec max (mm / day), 850 hpa temp anom (C))
                data_track  = np.ma.masked_all(13)
                data_track[0]   = time[time_i]
                data_track[1]   = lon_pres
                data_track[2]   = lat_pres
                data_track[3]   = pres_TC_day
                data_track[4]   = pres_TC_month
                data_track[5]   = vel_max
                data_track[6]   = rad_vel_max
                data_track[7]   = SST_TC_day
                data_track[8]   = SST_TC_month
                data_track[9]   = vor_850_TC_max
                data_track[10]  = TC_shear
                data_track[11]  = prec_max
                data_track[12]  = temp_anom_TC

                #Get the correct ID
                exec('TRACK_ID_'+str(track_i)+' = np.ma.vstack([TRACK_ID_'+str(track_i)+', data_track])')

                #Current RV maxima will not start a new track
                remove_index.append(index_min)
        PSL_min_lon     = np.delete(PSL_min_lon, remove_index)
        PSL_min_lat     = np.delete(PSL_min_lat, remove_index)
        temp_anom   = np.delete(temp_anom, remove_index)
        
        print(f"{counter_tooshort} active tracks (1) matched {len(remove_index)} valid")
        
        #-----------------------------------------------------------------------------------------
        # The remaining high RV areas will start new tracks
        
        for min_i in range(len(PSL_min_lat)):
            #Get the indices and coordinates for the PSL minimum
            lon_index   = PSL_min_lon[min_i]
            lat_index   = PSL_min_lat[min_i]
            lon_pres    = lon[lon_index]
            lat_pres    = lat[lat_index]

            #Get the 850 hPa TEMP anom (8x8)
            temp_anom_TC    = temp_anom[min_i]

            #Get the pressure and vorticity
            pres_TC_day = pres[time_i, lat_index, lon_index]
            pres_TC_month   = pres_month[lat_index, lon_index]

            #Get the maximum vorticity near the PSL minima
            vor_850_TC_max  = PressureRVMaxima(lon, lat, lon_index, lat_index, vor_850[time_i])

            #Find the radius of maximum wind speeds
            vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i])

            if rad_vel_max > 400:
                #Reduce region of search
                vel_max, rad_vel_max    = RadiusMaximumWind(lon, lat, lon_index, lat_index, U_10[time_i], reduced_radius = True)

            #Get the sea surface temperature
            SST_TC_day, SST_TC_month    = SeaSurfaceTemperatureTC(lon_ocn, lat_ocn, SST_day, SST_month, lon_pres, lat_pres)

            if SST_TC_day is np.ma.masked:
                #TC forms above land, discard RV maxima
                continue

            #Determine vertical wind shear
            TC_shear    = VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i], u_vel_250[time_i], v_vel_250[time_i])

            #Find the maximum precipitation
            prec_max    = MaximumPrecipitation(lon, lat, lon_index, lat_index, prec[time_i])

            #Data for track (time, lon (pres), lat (pres), pres (3 h), pres (month), vel_max, rad. vel_max, SST (day), SST (month), vor max (1/s), TC shear (m/s), prec max (mm / day), 850 hpa temp anom (C))
            data_track  = np.ma.masked_all(13)
            data_track[0]   = time[time_i]
            data_track[1]   = lon_pres
            data_track[2]   = lat_pres
            data_track[3]   = pres_TC_day
            data_track[4]   = pres_TC_month
            data_track[5]   = vel_max
            data_track[6]   = rad_vel_max
            data_track[7]   = SST_TC_day
            data_track[8]   = SST_TC_month
            data_track[9]   = vor_850_TC_max
            data_track[10]  = TC_shear
            data_track[11]  = prec_max
            data_track[12]  = temp_anom_TC

            #Initiate new track for new low
            exec('TRACK_ID_'+str(track_ID)+' = data_track')

            #Save all the ID's, some ID's will be removed (too short track etc.)
            track_ID_active.append(track_ID)
            
            #Raise track ID
            track_ID    += 1

        print(f"{len(PSL_min_lon)} new active tracks (1) created -> {len(track_ID_active)=}.")
        
        #-----------------------------------------------------------------------------------------
        # Now going through all active tracks, including the new ones and check if they are 
        # still active, i.e. if a new candidate point was found and it did not weaken twice in a row
        
        counter_notupdated = 0
        counter_tooshortregular = 0
        counter_tooshortland = 0
        counter_origin = 0
        counter_velmax = 0
        counter_shear = 0
        counter_SST = 0
        for track_i in track_ID_active:
            #Check the tracks whether they already finished
    
            #Get the last time stamp for each active track
            exec('track = TRACK_ID_'+str(track_i))

            if track.ndim == 1:
                #Only one time stamp
                time_end_track  = track[0]
                duration_track  = 3.0
        
            else:
                time_end_track  = track[-1, 0]
                duration_track  = len(track) * 3.0

            if time[time_i] == time_end_track:
                #Still an active track, keep active
                continue
            counter_notupdated += 1

            if duration_track >= 48.0:
                print(f"{track_i=} finished: {duration_track=}hr")
                #Storm lasts more than 2 days

                #Check whether the first 24 hours do not contain any land mask
                #Otherwise, remove this part of the track and conduct the checks
                for time_j in range(len(track) - 8):

                    if any(track[time_j:time_j + 8, 7].mask):
                        continue

                    #No adjustments needed, first 24 hour do not contain land mask
                    break

                if time_j != 0:
                    #Adjust the track and now conduct the final check
                    track   = track[time_j:]
                    exec('TRACK_ID_'+str(track_i)+' = track')

                    if len(track)  * 3.0 < 48:
                        #Duration of track is now less than 2 days
                        #Remove track (by setting start latitude at 90N)
                        print(f"{track_i=} too short after removing land time (set origin lat to 90N)")
                        counter_tooshortland += 1
                        track[0, 2] = 90


                vel_max_counter = 0
                shear_counter   = 0
                
                for time_j in range(len(track)):

                    if track[time_j, 5] >= 17.0:
                        #Check whether the velocity reaches at least 17 m/s
                        vel_max_counter += 1

                    else:
                        #Reset counter
                        vel_max_counter = 0

                    if vel_max_counter == 8:
                        #8 consecequtive time steps (24) of 17 m/s
                        break
                
                for time_j in range(len(track)):

                    if track[time_j, 10] <= 12.5:
                        #Check whether the shear is lower than 12.5 m/s
                        shear_counter += 1

                    else:
                        #Reset counter
                        shear_counter = 0

                    if shear_counter == 8:
                        #8 consecequtive time steps (24) of low shear
                        break
                        
                counter_origin += int(np.abs(track[0, 2]) > 30.0)
                counter_velmax += int(vel_max_counter != 8)
                counter_shear += int(shear_counter != 8)
                counter_SST += np.int64(track[0, 7] < 25.0) # 0 if masked
                
                if np.abs(track[0, 2]) <= 30.0 and vel_max_counter == 8 and shear_counter == 8 and track[0, 7] >= 25.0:
                    #Check whether TC originates between 30S and 30N
                    #Check whether max wind speed is sustained at 17 m/s for 24 hours
                    #Minimum sea surface temperature of 25C
                    track_ID_list.append(track_i)

                #Remove active track
                track_ID_active.remove(track_i)
            
            else:
                #Remove active track
                counter_tooshortregular += 1
                track_ID_active.remove(track_i)
#               print('line 1079: del TRACK_ID_'+str(track_i))
                exec('del TRACK_ID_'+str(track_i))

        print(f"{counter_notupdated} tracks finished -> {len(track_ID_active)=}\n  discard: {counter_tooshortregular+counter_tooshortland} short "
            + f"({counter_tooshortregular} reg. {counter_tooshortland=} land), {counter_origin} non trop., "
            + f"{counter_velmax} weak, {counter_shear} shear, {counter_SST} SST")
        print(f"  number of valid tracks: {len(track_ID_list)} (total)")
        #-----------------------------------------------------------------------------------------

# Now repeat for tracks that are still active at the end of the run
counter_tooshortregular = 0
counter_tooshortland = 0
counter_origin = 0
counter_velmax = 0
counter_shear = 0
counter_SST = 0
for track_i in track_ID_active:
    #Check the last tracks manually

    #Get the last time stamp for each active track
    exec('track = TRACK_ID_'+str(track_i))

    if track.ndim == 1:
        #Only one time stamp
        time_end_track  = track[0]
        duration_track  = 3.0

    else:
        time_end_track  = track[-1, 0]
        duration_track  = len(track) * 3.0

    if duration_track >= 48.0:
        #Storm lasts more than 2 days

        #Check whether the first 24 hours do not contain any land mask
        #Otherwise, remove this part of the track and conduct the checks
        for time_j in range(len(track) - 8):
            if any(track[time_j:time_j + 8, 7].mask):
                continue

            #No adjustments needed, first 24 hour do not contain land mask
            break

        if time_j != 0:
            #Adjust the track and now conduct the final check
            track   = track[time_j:]
            exec('TRACK_ID_'+str(track_i)+' = track')

            if len(track)  * 3.0 < 48:
                #Duration of track is now less than 2 days
                #Remove track from tracker
                counter_tooshortland += 1
                track[0, 2] = 90

        vel_max_counter = 0
        shear_counter   = 0
                
        for time_j in range(len(track)):

            if track[time_j, 5] >= 17.0:
                #Check whether the velocity reaches at least 17 m/s
                vel_max_counter += 1

            else:
                #Reset counter
                vel_max_counter = 0

            if vel_max_counter == 8:
                #8 consecequtive time steps (24) or 17 m/s
                break

        for time_j in range(len(track)):

            if track[time_j, 10] <= 12.5:
                #Check whether the shear is lower than 12.5 m/s
                shear_counter += 1

            else:
                #Reset counter
                shear_counter = 0

            if shear_counter == 8:
                #8 consecequtive time steps (24) of low shear
                break
                
        counter_origin += int(np.abs(track[0, 2]) > 30.0)
        counter_velmax += int(vel_max_counter != 8)
        counter_shear += int(shear_counter != 8)
        counter_SST += np.int64(track[0, 7] < 25.0)
        
        if np.abs(track[0, 2]) <= 30.0 and vel_max_counter == 8 and shear_counter == 8 and track[0, 7] >= 25.0:
            #Check whether TC originates between 30S and 30N
            #Check whether max wind speed is sustained at 17 m/s for 24 hours
            #Minimum sea surface temperature of 25C
            print(f"Track {track_i} is complete and valid.")
            track_ID_list.append(track_i)
            
    else:
        counter_tooshortregular += 1
print(f"Tracks still active at the end ({len(track_ID_active)}):\n  {counter_tooshortregular+counter_tooshortland} short "
    + f"({counter_tooshortregular} reg. {counter_tooshortland=} land), {counter_origin} non trop., "
    + f"{counter_velmax} weak, {counter_shear} shear, {counter_SST} SST")       
print(f"Total number of valid stored tracks: {len(track_ID_list)}")
#-----------------------------------------------------------------------------------------

#Determine the maximum track time
time_max    = 0

for track_i in track_ID_list:

    exec('track = TRACK_ID_'+str(track_i))

    if len(track) > time_max:
        #Update maximum time
        time_max    = len(track)

#-----------------------------------------------------------------------------------------

#Raise time max by 1 (easy to remove masked objects for post-processing)
time_max    += 1

#Save all the active tracks
track_all   = np.ma.masked_all((len(track_ID_list), time_max, 13))

for track_i in range(len(track_ID_list)):
    #Save data to general array
    exec('track = TRACK_ID_'+str(track_ID_list[track_i]))
    track_all[track_i, :len(track)] = track
    
track_all[:,:,7:9][track_all[:,:,7:9]==100] = np.nan # missing SST/SSTmon values

#-----------------------------------------------------------------------------------------
#Save the coordinates for each track
HEAT_data = netcdf.Dataset(f'{directory}/{outfile}', 'w')

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
print(f"Results saved to {directory}/{outfile}")
