#!/usr/bin/python3.6

#Program determines all the relative vorticity maxima

import cftime
import numpy as np
import datetime
import time
import glob, os
import sys
import math
import netCDF4 as netcdf
import multiprocessing as mp

dry_run			= False # set False to do computing

experiment = sys.argv[1] # run with e.g. python RVmax_finder.py rcp 1
run_ensemble = int(sys.argv[2])

#-----------------------------------------------------------------------------------------------------------------
# uncomment for RCP8.5 hres data
if experiment == 'ref':
	run_year_start	= 2002 # [2002, 2092]
	experiment_name = 'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)
	directory_data  = f'/projects/0/prace_imau/prace_2013081679/cesm1_0_4/f02_t12/{experiment_name}/OUTPUT/atm/hist/3h/'
	stream		  	= 'h1'
	directory		= f'/home/jasperdj/files_rene/RCP.started_{run_year_start}.{run_ensemble:03d}/RV_Max_seeds/'
	gridfile		= '/home/jasperdj/files_rene/Atmosphere_0_25_DX_DY_AREA.nc'
	NOUT			= 8 # number of time steps per output file (one day)
elif experiment == 'rcp':
	run_year_start	= 2092 # [2002, 2092]
	experiment_name = 'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)
	directory_data  = f'/projects/0/prace_imau/prace_2013081679/cesm1_0_4/f02_t12/{experiment_name}/OUTPUT/atm/hist/3h/'
	stream		  	= 'h1'
	directory		= f'/home/jasperdj/files_rene/RCP.started_{run_year_start}.{run_ensemble:03d}/RV_Max_seeds/'
	gridfile		= '/home/jasperdj/files_rene/Atmosphere_0_25_DX_DY_AREA.nc'
	NOUT			= 8 # number of time steps per output file (one day)
elif experiment == 'sai':
	run_year_start	= 2092 # [2092]
	experiment_name = 'hres_b.e10.B2000_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)
	directory_data  = '/projects/0/nwo2021025/archive/'+experiment_name+'/atm/hist/'
	stream		  	= 'h5'
	directory		= f'/home/jasperdj/files_rene/SAI.started_{run_year_start}.{run_ensemble:03d}/RV_Max_seeds/'
	gridfile		= '/home/jasperdj/files_rene/Atmosphere_0_25_DX_DY_AREA.nc'
	NOUT			= 8 # number of time steps per output file (one day)
else:
	raise ValueError(f'Wrong input argument: {experiment=}')

# # uncomment for SAI mres (0.5degree) data
# experiment_name = 'mres_b.e10.B2000_CAM5.f05_t12.001'
# directory_data	= '/projects/0/nwo2021025/archive/'+experiment_name+'/atm/hist/'
# stream			= 'h4'
# stream3D		= 'h3' # only used for mres to calculate U250, V250 and T850
# directory		= '/home/jasperdj/files_rene/SAI.mres/RV_Max/'
# gridfile		= '/home/jasperdj/files_rene/Atmosphere_0_5_DX_DY_AREA.nc'
# NOUT			= 4 # number of time steps per output file (one day)
#-----------------------------------------------------------------------------------------------------------------

# input data
files = glob.glob(directory_data+experiment_name+'.cam2.'+stream+'.*.nc')
files.sort() #Sort the files on date
# files = files[6:372]

LATMIN = -60
LATMAX = 60 
RADIUSEARTH=6371000.0
RVTHRESHOLD=6.0 * 10**(-5.0) # minimum for RV@850hPa
#RVDIFFTHRESHOLD=6.0 * 10**(-5.0) # minimum for RV@850hPa - RV@250hPa
#CORETEMPTHRESHOLD=0.0 # minimum temperature anomaly core w.r.t. 8x8 degree area
FIELDSIZE=8 # (degrees N x E) size of field to calculate reference temperature
#U10THRESHOLD=10 # (m/s) minimum 10m wind speed
#U10THRESHOLD_MAXRADIUS=100 # (km) U10 should be 10m/s within 100km
RVDISTMIN=250 # (km) min. distance between maxima (if less, only the one with lower PSL counts)


def chunk_data(files):
	"""divide data into daily chunks
	note: date is determined by center of time_bnds, not time!
	
	result: chunks = [(date, fileids, t0, t1), ...]
		such that netcdf.MFDataset(files[fileids])['time'][t0:t1] gives
		the timestamps corresponding to date
	"""
	n = 0 # counter for time step
	with netcdf.Dataset(files[0],'r') as fh:
# 		time = fh.variables['time'][:]
		time = fh['time_bnds'][:].mean(axis=-1)
# 		if 'time: mean' in # getattr(fh['U850'], 'cell_methods', ''):
# 			time = time - 0.5*dt
		dates = [cftime.num2date(t, fh['time'].units, fh['time'].calendar).strftime("%Y-%m-%d") for t in time]
		chunks = [[dates[0],[0],0,0]]
	for fid in range(len(files)):
		with netcdf.Dataset(files[fid],'r') as fh:
# 			time = fh.variables['time'][:]
			time = fh['time_bnds'][:].mean(axis=-1)
# 			if 'time: mean' in getattr(fh['U850'], 'cell_methods', ''):
# 				time = time - 0.5*dt
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


# # old version
# def chunk_data(files):
# 	"""divide data into daily chunks"""
# 	n = 0
# 	with netcdf.Dataset(files[0],'r') as fh:
# 		Time = fh.variables['time']
# 		CTime = fh.variables['time_bnds'][:].mean(axis=1)
# 		dates = [cftime.num2date(T, Time.units, Time.calendar).strftime("%Y%m%d") for T in CTime]
# 		chunks = [[dates[0],[0],0,0]]
# 	for fid in range(len(files)):
# 		with netcdf.Dataset(files[fid],'r') as fh:
# 			CTime = fh.variables['time_bnds'][:].mean(axis=1)
# 			dates = [cftime.num2date(T, Time.units, Time.calendar).strftime("%Y%m%d") for T in CTime]
# 			for i,newdate in enumerate(dates):
# 				(olddate, fids0, i0, n0) = chunks[-1]
# 				if newdate != olddate:
# 					if fids0[-1] != fid and i!=0:
# 						chunks[-1][1].append(fid)
# 					chunks[-1][3] = i0 + n - n0
# 					chunks.append([newdate,[fid],i,n])
# 				n += 1
# 	(olddate, fids0, i0, n0) = chunks[-1]
# 	chunks[-1][3] = i0 + n - n0
# 	return chunks


# def ReadinDataFixForMRES(filenames, t1, t2):
# 	"""The mres h4 file does not contain U250, V250 and T850, so using 3D data"""
# 	with netcdf.MFDataset(filenames, 'r', aggdim='time') as fh:
# 		time	   = fh.variables['time']				   #Time
# 		time	   = cftime.num2date(time[t1:t2], time.units, calendar=time.calendar)
# 		lon			= fh.variables['lon'][:]					# Longitude
# 		lat			= fh.variables['lat'][i1-1:i2+1]			# Latitude
# 		lat_weight	= fh.variables['gw'][i1:i2]				    # Latitude weight
# 		pres		= fh.variables['PSL'][t1:t2, i1:i2] / 100.0	# Sea level pressure (hPa)
# 		U_10		= fh.variables['U10'][t1:t2, i1:i2]		    # 10-meter wind speed (m/s)
# 		u_vel_850	= fh.variables['U850'][t1:t2, i1-1:i2+1]    # Zonal velocity at 850 hPa (m/s)
# 		v_vel_850	= fh.variables['V850'][t1:t2, i1-1:i2+1]    # Zonal velocity at 850 hPa (m/s)
# 	files3D = glob.glob(directory_data+experiment_name+'.cam2.'+stream3D+'.*.nc')
# 	files3D.sort() #Sort the files on date
# 	fidi = files3D.index(filenames[0].replace(f".{stream}.", f".{stream3D}."))
# 	for fid in range(fidi,len(files3D)):
# 		with netcdf.Dataset(files3D[fid], 'r') as fh:
# 			time3D = fh.variables['time']
# 			time3D = cftime.num2date(time3D[:],time3D.units,calendar=time3D.calendar)
# 		if time[0] in time3D:
# 			fid0 = fid
# 		if time[-1] in time3D:
# 			fid1 = fid
# 			break
# 	with netcdf.MFDataset(files3D[fid0:fid1+1], 'r', aggdim='time') as fh:
# 		time3D = fh.variables['time']
# 		time3D = cftime.num2date(time3D[:],time3D.units,calendar=time3D.calendar)
# 		ti, tf = time3D.searchsorted([time[0], time[-1]])
# 		assert all(time3D[ti:tf+1]==time), (
# 			f"Could not obtain steps {time=} from {time3D=} in {stream3D} output.")
# 		lev = fh.variables['lev'][:]
# 		assert max(lev) < 1200, f"cannot assume lev has units hPa (max: {max(lev)=})" # assume hPa unit
# 		l250 = abs(lev-250).argmin()
# 		l850 = abs(lev-850).argmin()
# 		print(f"Warning @ {[os.path.basename(f) for f in filenames]}\n  reading U, V at lev={lev[l250]:.1f}hPa, "
# 			+ f"T at lev={lev[l850]:.1f}hPa instead of true pressure levels.")
# 		temp = fh.variables['T'][ti:tf+1, l850, i1:i2]
# 		u_vel_250 = fh.variables['U'][ti:tf+1, l250, i1-1:i2+1]
# 		v_vel_250 = fh.variables['V'][ti:tf+1, l250, i1-1:i2+1]
# 		return time, lon, lat, lat_weight, pres, temp, U_10, u_vel_850, u_vel_250, v_vel_850, v_vel_250


def ReadinData(filenames, x_diff, y_diff, t1, t2):
	# if experiment_name == 'mres_b.e10.B2000_CAM5.f05_t12.001':
# 		(time, lon, lat, lat_weight, pres, temp, U_10, u_vel_850,
# 		u_vel_250, v_vel_850, v_vel_250) = ReadinDataFixForMRES(filenames, t1, t2)
# 	else:
	with netcdf.MFDataset(filenames, 'r', aggdim='time') as fh:
		time	   = fh.variables['time'][t1:t2]				   #Time (end of interval)
		if 'time: mean' in getattr(fh['U850'], 'cell_methods', ''):
			time = fh.variables['time_bnds'][:].mean(axis=-1) # only center data is time averaged
		time	   = cftime.num2date(time, fh['time'].units, fh['time'].calendar)
		date	   = fh.variables['time_bnds'][t1:t2].mean(axis=-1) # date always uses centered time
		date 	   = np.array([cftime.num2date(d, fh['time'].units, fh['time'].calendar).strftime('%Y-%m-%d') for d in date])
		lon		   = fh.variables['lon'][:]					   #Longitude
		lat		   = fh.variables['lat'][i1-1:i2+1]				 #Latitude
		lat_weight = fh.variables['gw'][i1:i2]				 #Latitude weight
		pres	   = fh.variables['PSL'][t1:t2, i1:i2] / 100.0	 #Sea level pressure (hPa)
		U_10	   = fh.variables['U10'][t1:t2, i1:i2]		 #10-meter wind speed (m/s)
		u_vel_850  = fh.variables['U850'][t1:t2, i1-1:i2+1]	 #Zonal velocity at 850 hPa (m/s)
		v_vel_850  = fh.variables['V850'][t1:t2, i1-1:i2+1]	 #Zonal velocity at 850 hPa (m/s)
		temp	   = fh.variables['T850'][t1:t2, i1:i2]		 #Temperature 850 hPa
		u_vel_250  = fh.variables['U250'][t1:t2, i1-1:i2+1]	 #Zonal velocity at 250 hPa (m/s)
		v_vel_250  = fh.variables['V250'][t1:t2, i1-1:i2+1]	 #Zonal velocity at 250 hPa (m/s)

	#Add periodic boundaries
	lon_2, u_vel_850	= PeriodicBoundaries3D(time, lon, lat, u_vel_850)
	lon_2, v_vel_850	= PeriodicBoundaries3D(time, lon, lat, v_vel_850)
	lon_2, u_vel_250	= PeriodicBoundaries3D(time, lon, lat, u_vel_250)
	lon, v_vel_250		= PeriodicBoundaries3D(time, lon, lat, v_vel_250)

	#Determine the vorticity at 850 and 250 hPa
	vor_850			= np.ma.masked_all((len(time), len(lat) - 2, len(lon) - 2))
	vor_250			= np.ma.masked_all((len(time), len(lat) - 2, len(lon) - 2))

	for time_i in range(len(time)):
		#Vorticity for the two level
		vor_850[time_i] = RelativeVorticity(u_vel_850[time_i], v_vel_850[time_i], x_diff, y_diff)
		vor_250[time_i] = RelativeVorticity(u_vel_250[time_i], v_vel_250[time_i], x_diff, y_diff)
	
	#Get the correct lon/lat dimensions
	lon, lat		= lon[1:-1], lat[1:-1]

	return time, date, lon, lat, lat_weight, pres, U_10, temp, vor_850, vor_250


def PeriodicBoundaries2D(lon, lat, field, lon_grids = 1):
	"""Add periodic zonal boundaries for 2D field"""

	#Empty field with additional zonal boundaries
	lon_2			= np.zeros(len(lon) + lon_grids * 2)
	field_2			= np.ma.masked_all((len(lat), len(lon_2)))
	
	#Get the left boundary, which is the right boundary of the original field
	lon_2[:lon_grids]	= lon[-lon_grids:] - 360.0
	field_2[:, :lon_grids]	= field[:, -lon_grids:]

	#Same for the right boundary
	lon_2[-lon_grids:]	= lon[:lon_grids] + 360.0
	field_2[:, -lon_grids:] = field[:, :lon_grids]

	#And the complete field
	lon_2[lon_grids:-lon_grids]		= lon
	field_2[:, lon_grids:-lon_grids]	= field

	return lon_2, field_2	


def PeriodicBoundaries3D(time, lon, lat, field, lon_grids = 1):
	"""Add periodic zonal boundaries for 3D field"""

	#Empty field with additional zonal boundaries
	lon_2				= np.zeros(len(lon) + lon_grids * 2)
	field_2				= np.ma.masked_all((len(time), len(lat), len(lon_2)))
	
	#Get the left boundary, which is the right boundary of the original field
	lon_2[:lon_grids]		= lon[-lon_grids:] - 360.0
	field_2[:, :, :lon_grids]	= field[:, :, -lon_grids:]

	#Same for the right boundary
	lon_2[-lon_grids:]		= lon[:lon_grids] + 360.0
	field_2[:, :, -lon_grids:]	= field[:, :, :lon_grids]

	#And the complete field
	lon_2[lon_grids:-lon_grids]		= lon
	field_2[:, :, lon_grids:-lon_grids]		= field

	return lon_2, field_2	


def RelativeVorticity(u_vel, v_vel, x_diff, y_diff):
	"""Determines the relative vorticity of the field"""

	#Take the meridional difference of the zonal wind
	u_vel_diff	= u_vel[2:] - u_vel[:-2]
	u_vel_diff	= u_vel_diff[:, 1:-1]

	#Take the zonal difference of the meridional wind
	v_vel_diff	= v_vel[:, 2:] - v_vel[:, :-2]
	v_vel_diff	= v_vel_diff[1:-1]

	#Determine the relative vorticity
	vorticity	= (v_vel_diff / x_diff) - (u_vel_diff / y_diff)

	return vorticity


def Distance(lon_1, lat_1, lon_2, lat_2):
	"""Returns distance (m) of two points located at the globe coordinates need input in degrees"""

	#Convert to radians
	lon_1, lat_1, lon_2, lat_2 = np.radians([lon_1, lat_1, lon_2, lat_2]) 

	#Haversine formula 
	d_lon	= lon_2 - lon_1 
	d_lat	= lat_2 - lat_1 
	a	= math.sin(d_lat/2.0)**2 + math.cos(lat_1) * math.cos(lat_2) * math.sin(d_lon/2.0)**2
	c	= 2.0 * math.asin(np.sqrt(a)) 
	
	return c * RADIUSEARTH #Distance between two points in meter


def main(task):
	checkdate, filelist, t0, t1 = task
	try:
		filenames = [files[fid] for fid in filelist]
	except IndexError:
		print(f"IndexError: task {task} failed. {len(files)=}")
		return
	#print(f"[{datetime.datetime.now()}] File {filelist}, steps {t0} - {t1}", flush=True)
	
	#Counters
	hour_counter		= 0
	RV_max_max_day		= 0

	fh = netcdf.Dataset(gridfile, 'r')

	#Writing data to correct variable
	lon_grid	= fh.variables['lon'][:]	#Longitude
	lat_grid	= fh.variables['lat'][i1-1:i2+1]  #Latitude	  
	grid_x		= fh.variables['DX'][i1-1:i2+1]	  #Zonal length of grid cell  
	grid_y		= fh.variables['DY'][i1-1:i2+1]	  #Meridional length of grid cell	  

	fh.close()

	lon_2, grid_x	= PeriodicBoundaries2D(lon_grid, lat_grid, grid_x)
	lon, grid_y = PeriodicBoundaries2D(lon_grid, lat_grid, grid_y)

	#Generate for each grid cell the difference length
	x_diff		= 0.5 * grid_x[:, 2:] + 0.5 * grid_x[:, :-2] + grid_x[:, 1:-1]
	y_diff		= 0.5 * grid_y[2:] + 0.5 * grid_y[:-2] + grid_y[1:-1]
	x_diff	= x_diff[1:-1]
	y_diff	= y_diff[:, 1:-1]
	
	#For each day (8x3 hour) determine the pressure lows
	time_RV_max		= np.ma.masked_all(NOUT)
	RV_max_coor_lon		= np.ma.masked_all((NOUT, 800))
	RV_max_coor_lat		= np.ma.masked_all((NOUT, 800))
	temp_anom_all		= np.ma.masked_all((NOUT, 800))
	time, date, lon, lat, lat_weight, pres, U_10, temp, vor_850, vor_250	= ReadinData(filenames, x_diff, y_diff, t0, t1)
	assert (lon_grid.shape == lon.shape) and (lat_grid[1:-1].shape == lat.shape), (
		f'Grid in {gridfile} does not match data.')
		
	# set RV850 to 0 outside the tropics for TC seeds
	vor_850 = np.where(np.abs(lat[np.newaxis,:,np.newaxis]) <= 30, vor_850, 0)

	#Last elements can not be determined
	x_diff	= x_diff[1:-1]
	y_diff	= y_diff[:, 1:-1]
   
	for time_i in range(len(time)):
		date_str = date[time_i] #time[time_i].strftime('%Y-%m-%d')
		if os.path.exists(directory+'RV_Max_Coordinates_'+date_str+'.nc'):
			hour_counter	= 0
			#print('skip creating RV_Max_Coordinates_'+date_str+'.nc (already exists)', flush=True)
			continue
		else:
			print(f"[{datetime.datetime.now()}] MFFile ids: {filelist}, step: {time_i+t0:03d} -> {date_str} ({t0:03d}-{t1:03d})", flush=True)
		assert date_str == checkdate, f"{date_str=} not equal to {checkdate=}"
		
		#Loop over each 3-hourly field
		lat_index_EQ	= np.where(lat >= 0.0)[0][0]

		#Get all the indices where the (absolute) vorticity is larger than 6 * 10^-5 
		index_NH	= np.where(vor_850[time_i, lat_index_EQ:] >= RVTHRESHOLD)
		lat_index_NH	= index_NH[0] + lat_index_EQ
		lon_index_NH	= index_NH[1]
		index_SH	= np.where(-vor_850[time_i, :lat_index_EQ] >= RVTHRESHOLD)
		lat_index_SH	= index_SH[0]
		lon_index_SH	= index_SH[1]

		#Save all the indices in 1 array (SH + NH)
		lon_index	= np.insert(lon_index_SH, len(lon_index_SH), lon_index_NH)
		lat_index	= np.insert(lat_index_SH, len(lat_index_SH), lat_index_NH)

		#Save all the indices in 1 array (SH + NH)
		lon_index	= np.insert(lon_index_SH, len(lon_index_SH), lon_index_NH)
		lat_index	= np.insert(lat_index_SH, len(lat_index_SH), lat_index_NH)

		#-----------------------------------------------------------------------------------------
		#Empty arrays to save the RV maxima
		RV_max_index_lon	= []
		RV_max_index_lat	= []
		temp_anom		= []

		for index_i in range(len(lat_index)):
			#Loop over each point with a negative sea-level pressure
			lon_i		= lon_index[index_i]
			lat_i		= lat_index[index_i]

			#Check vorticity criterium
			vor_850_TC	= vor_850[time_i, lat_i, lon_i] * np.sign(lat[lat_i])
# 			vor_250_TC	= vor_250[time_i, lat_i, lon_i] * np.sign(lat[lat_i])

			if vor_850_TC < RVTHRESHOLD:
				#Relative vorticity is too low
				continue
			
			assert (np.abs(lat[lat_i])<=30), f"{checkdate=}, {time_i=}, lat={lat[lat_i]}, RV={vor_850[time_i,lat_i,lon_i]}" # <<<<<<<<<<<<<<<

# 			if (vor_850_TC - vor_250_TC) < RVDIFFTHRESHOLD:
# 				#No evidence for warm core
# 				continue

			#Get the horizontal velocity (8x8) around the eye
			lon_east	= lon_i + round(FIELDSIZE/(2*dlon)) + 1
			lon_west	= lon_i - round(FIELDSIZE/(2*dlon))
			lat_north	= lat_i + round(FIELDSIZE/(2*dlat)) + 1
			lat_south	= lat_i - round(FIELDSIZE/(2*dlat))

			if lat_south < 0:
				lat_south = 0

			#Get the dimensions for the lat field
			lat_field	= lat[lat_south:lat_north]
			weight		= lat_weight[lat_south:lat_north]
			weight		= weight / np.sum(weight)
			lon_field	= np.zeros(lon_east - lon_west)
			temp_field	= np.ma.masked_all((len(lat_field), len(lon_field)))
			U_10_field	= np.ma.masked_all((len(lat_field), len(lon_field)))

			if lon_east > len(lon):
				#Eastern part is on eastern hemisphere
				lon_1				= lon[lon_west:]
				lon_2				= lon[:lon_east - len(lon)] + 360.0
				lon_field[:len(lon_1)]		= lon_1
				lon_field[len(lon_1):]		= lon_2
				temp_field[:, :len(lon_1)]	= temp[time_i, lat_south:lat_north, lon_west:]
				temp_field[:, len(lon_1):]	= temp[time_i, lat_south:lat_north, :lon_east - len(lon)]
				U_10_field[:, :len(lon_1)]	= U_10[time_i, lat_south:lat_north, lon_west:]
				U_10_field[:, len(lon_1):]	= U_10[time_i, lat_south:lat_north, :lon_east - len(lon)]

			elif lon_west < 0:
				#Western part is on western hemisphere
				lon_1				= lon[lon_west:] - 360.0
				lon_2				= lon[:lon_east]
				lon_field[:len(lon_1)]		= lon_1
				lon_field[len(lon_1):]		= lon_2
				temp_field[:, :len(lon_1)]	= temp[time_i, lat_south:lat_north, lon_west:]
				temp_field[:, len(lon_1):]	= temp[time_i, lat_south:lat_north, :lon_east]
				U_10_field[:, :len(lon_1)]	= U_10[time_i, lat_south:lat_north, lon_west:]
				U_10_field[:, len(lon_1):]	= U_10[time_i, lat_south:lat_north, :lon_east]

			else:
				#Normal field can be retained
				lon_field	= lon[lon_west:lon_east]
				temp_field	= temp[time_i, lat_south:lat_north, lon_west:lon_east]
				U_10_field	= U_10[time_i, lat_south:lat_north, lon_west:lon_east]

			#Determine the spatial average over the field
			temp_mean	= np.mean(temp_field, axis = 1)
			temp_mean	= np.sum(temp_mean * weight)

			#Get the temperature anomaly of the given point
			lon_index_2 = (np.abs(lon_field - lon[lon_i])).argmin()
			lat_index_2 = (np.abs(lat_field - lat[lat_i])).argmin()
			temp_anom_TC	= temp_field[lat_index_2, lon_index_2] - temp_mean

# 			if temp_anom_TC < CORETEMPTHRESHOLD:
# 				#No warm core
# 				continue
# 
# 			#Check whether the 10-m wind speed exceeds threshold 
# 			U_10_treshold = False
# 
# 			for lat_j in range(len(lat_field)):
# 				for lon_j in range(len(lon_field)):
# 					#Check whether the 10-m wind speeds exceeds 10 m/s
# 					if U_10_field[lat_j, lon_j] < U10THRESHOLD:
# 						continue
# 
# 					#Check the 10-m wind speed in a 100-km search radius
# 					distance	= Distance(lon[lon_i], lat[lat_i], lon_field[lon_j], lat_field[lat_j]) / 1000.0
# 
# 					if distance <= U10THRESHOLD_MAXRADIUS:
# 						#10-m wind speed exceeds treshold within 100 km
# 						U_10_treshold = True
# 						break
# 
# 				if U_10_treshold:
# 					#Threshold is reached, stop search
# 					break
# 
# 			if U_10_treshold == False:
# 				#The 10-m wind speed treshold is not valid
# 				continue

			#Save the coordinates and the anomalies
			RV_max_index_lon.append(lon_i)
			RV_max_index_lat.append(lat_i)
			temp_anom.append(temp_anom_TC)

		#-----------------------------------------------------------------------------------------
		distance_RV = True
		RV_index_start	= 0

		while distance_RV:
			#Check RV maxima in the neighbourhood
			RV_dlat_min = np.ceil(1000*RVDISTMIN*180/(np.pi*RADIUSEARTH))
			for RV_i in range(RV_index_start, len(RV_max_index_lon)):
				for RV_j in range(RV_i + 1, len(RV_max_index_lon)):
					#Get the coordinates
					lon_RV_1	= lon[RV_max_index_lon[RV_i]]
					lat_RV_1	= lat[RV_max_index_lat[RV_i]]
					lon_RV_2	= lon[RV_max_index_lon[RV_j]]
					lat_RV_2	= lat[RV_max_index_lat[RV_j]]

					if np.abs(lat_RV_1 - lat_RV_2) > RV_dlat_min:
						#Points are too far apart
						continue

					#Determine the distance (km) between the two points
					distance	= Distance(lon_RV_1, lat_RV_1, lon_RV_2, lat_RV_2) / 1000.0

					if distance < RVDISTMIN:
						#Distance between pressure minima must be at least 250 km
						#Select the one with the lowest sea-level pressure
						pres_1		= pres[time_i, RV_max_index_lat[RV_i], RV_max_index_lon[RV_i]]
						pres_2		= pres[time_i, RV_max_index_lat[RV_j], RV_max_index_lon[RV_j]]

						#Retain the one with the highest RV
						index_min	= np.argmin([pres_1, pres_2])

						if index_min == 0:
							#First point has lowest PSL, remove second
							RV_max_index_lon	= np.delete(RV_max_index_lon, RV_j)
							RV_max_index_lat	= np.delete(RV_max_index_lat, RV_j)
							temp_anom		= np.delete(temp_anom, RV_j)

						else:
							#Second point has lowest PSL, remove first
							RV_max_index_lon	= np.delete(RV_max_index_lon, RV_i)
							RV_max_index_lat	= np.delete(RV_max_index_lat, RV_i)
							temp_anom		= np.delete(temp_anom, RV_i)

						#Break the j-loop
						RV_j = -1
						break

				#Check whether the j-loop was terminated
				if RV_j == -1:
					#Break i-loop and start over with the reduced list
					RV_j = 0
					break

				else:
					#j-loop is completed, raise starting index by 1
					RV_index_start += 1

				if RV_i == len(RV_max_index_lon) - 1:
					#i-loop is finished, stop the while
					distance_RV = False


		#-----------------------------------------------------------------------------------------

		#Save the indices (note that the index has a base starting at 60S)
		time_RV_max[hour_counter]				= cftime.date2num(time[time_i], 'Days since 0001-01-01 00:00:00 UTC', 'noleap')
		RV_max_coor_lon[hour_counter, :len(RV_max_index_lon)]	= RV_max_index_lon
		RV_max_coor_lat[hour_counter, :len(RV_max_index_lat)]	= RV_max_index_lat
		temp_anom_all[hour_counter, :len(temp_anom)]		= temp_anom

		if len(RV_max_index_lon) > RV_max_max_day:
			#New hourly maximum of pressure minima
			RV_max_max_day = len(RV_max_index_lon)

		#Update hour counter
		hour_counter += 1

		if hour_counter == NOUT:
			#Remove masked arrays
			RV_max_coor_lon		= RV_max_coor_lon[:, :RV_max_max_day + 1]
			RV_max_coor_lat		= RV_max_coor_lat[:, :RV_max_max_day + 1]
			temp_anom_all		= temp_anom_all[:, :RV_max_max_day + 1]

			# write data
			HEAT_data = netcdf.Dataset(directory+'RV_Max_Coordinates_'+date_str+'.nc', 'w')
			HEAT_data.createDimension('time', len(time_RV_max))
			HEAT_data.createDimension('number_lows', len(RV_max_coor_lon[0]))
			HEAT_data.createVariable('time', float, ('time'), zlib=True)
			HEAT_data.createVariable('number_lows', float, ('number_lows'), zlib=True)
			HEAT_data.createVariable('lon_index', float, ('time', 'number_lows'), zlib=True)
			HEAT_data.createVariable('lat_index', float, ('time', 'number_lows'), zlib=True)
			HEAT_data.createVariable('TEMP', float, ('time', 'number_lows'), zlib=True)
			HEAT_data.variables['lon_index'].longname	= 'Longitude index of low'
			HEAT_data.variables['lat_index'].longname	= 'Latitude index of low (0 = 60.2S)'
			HEAT_data.variables['TEMP'].longname		= 'Temperature anomaly w.r.t. 8x8 mean'
			HEAT_data.variables['time'].units		= 'Days since 0001-01-01 00:00:00 UTC'
			HEAT_data.variables['time'].calendar	= 'noleap'
			HEAT_data.variables['TEMP'].units		= 'deg C'	
			HEAT_data.variables['time'][:]			= time_RV_max
			HEAT_data.variables['number_lows'][:]	= np.arange(len(RV_max_coor_lon[0])) + 1
			HEAT_data.variables['lon_index'][:]		= RV_max_coor_lon
			HEAT_data.variables['lat_index'][:]		= RV_max_coor_lat
			HEAT_data.variables['TEMP'][:]			= temp_anom_all
			print("Created {}".format(HEAT_data.filepath()), flush=True)
			HEAT_data.close()
			

			#Reset the counters
			hour_counter	= 0
			RV_max_max_day	= 0
		
			#Make empty array for the following day
			time_RV_max		= np.ma.masked_all(NOUT)
			RV_max_coor_lon = np.ma.masked_all((NOUT, 800))
			RV_max_coor_lat = np.ma.masked_all((NOUT, 800))
			temp_anom_all	= np.ma.masked_all((NOUT, 800))
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

# only use data from 60S to 60N
try:
	with netcdf.Dataset(files[0], 'r') as fh:
		lats0 = fh['lat'][:]
		lons0 = fh['lon'][:]
		timei = fh['time'][0]
		dt = fh['time'][1] - fh['time'][0]
		assert np.all(lats0[1:] > lats0[:-1]), "latitude does not increase"
		i1 = lats0.searchsorted(LATMIN) - 1 # lat. index of 60S
		i2 = lats0.searchsorted(LATMAX) + 1 # lat. index of 60N
		dlat = np.mean(np.diff(lats0)) # meridional grid cell size [degrees N]
		dlon = np.mean(np.diff(lons0)) # zonal grid cell size [degrees E]
	with netcdf.Dataset(files[-1], 'r') as fh:
		try: # get last timestamp 
			timef = fh['time'][-1]
		except IndexError: # read one before last file if last file is empty
			with netcdf.Dataset(files[-2],'r') as fh2:
				timef = fh2['time'][-1]
		total_num_time_steps = (timef - timei)//dt + 1
except IndexError:
	print("could not find any files with syntax: "+directory_data+experiment_name+".cam2."+stream+".*.nc")
	raise


if __name__ == '__main__':
	print('output directory:', directory)
	print('First file: '+files[0])
	print('Last file: '+files[-1])
	print(f'{len(files)=}')
	print(f"using lat slice: lat[{i1}:{i2}]")
	print(f"{dlat=:.4f}, {dlon=:.4f}")
	
	# divide into manageable chunks of data
	estimated_num_total_tasks = total_num_time_steps // NOUT
	num_files = len(glob.glob(directory+'*.nc'))
	estimated_num_unfinished_tasks = estimated_num_total_tasks - num_files
	print(f"{estimated_num_total_tasks=}")
	print(f"{estimated_num_unfinished_tasks=}")

	# create output paths
	if os.path.exists(directory):
		print(f"Writing to {directory}")
	else:
		print(f"Creating {directory}")
		if not dry_run:
			os.makedirs(directory)

	# execute main function in parallel
	ncores = int(os.environ.get('SLURM_JOB_CPUS_PER_NODE',1))
	print(f"Starting processes on {ncores=}")
	if dry_run:
		sys.exit()
	mp.set_start_method('spawn')
	print("making chunks...")
	tasks = chunk_data(files)
	print("starting computations...")
	with mp.Pool(ncores) as pool:
		pool.map(main, tasks)

