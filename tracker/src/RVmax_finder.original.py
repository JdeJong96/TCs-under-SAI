#Program determines all the relative vorticity maxima

from matplotlib import pylab
from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
#directory 		= '/home/rvwesten/Tropical_Cyclones/Data/RCP/'
directory 		= '/home/jasperdj/files_rene/'
run_year_start		= 2002
run_ensemble		= 1 
directory_data 		= '/projects/0/prace_imau/prace_2013081679/cesm1_0_4/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'/OUTPUT/atm/hist/3h/'

def ReadinData(filename, x_diff, y_diff):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		#Time
	lon		= fh.variables['lon'][:]		#Longitude
	lat		= fh.variables['lat'][126:642]		#Latitude
	lat_weight	= fh.variables['gw'][127:641]		#Latitude weight
	pres    	= fh.variables['PSL'][:, 127:641] / 100.0	#Sea level pressure (hPa)
	U_10	  	= fh.variables['U10'][:, 127:641]	#10-meter wind speed (m/s)
	temp	  	= fh.variables['T850'][:, 127:641]	#Temperature 850 hPa
	u_vel_850    	= fh.variables['U850'][:, 126:642]	#Zonal velocity at 850 hPa (m/s)
	v_vel_850    	= fh.variables['V850'][:, 126:642]	#Zonal velocity at 850 hPa (m/s)
	u_vel_250    	= fh.variables['U250'][:, 126:642]	#Zonal velocity at 250 hPa (m/s)
	v_vel_250    	= fh.variables['V250'][:, 126:642]	#Zonal velocity at 250 hPa (m/s)

	fh.close()

	#Set the time to toordinal form
	date  		= filename[-19:-9]	
	year  		= int(date[0:4])
	month 		= int(date[5:7])
	day		= int(date[8:10])
	
	#Get the difference (3 hours in days)
	time_diff	= time - time[0]
	time		= datetime.datetime(year, month, day).toordinal() + time_diff + time_diff[1] / 2.0

	#Add periodic boundaries
	lon_2, u_vel_850	= PeriodicBoundaries3D(time, lon, lat, u_vel_850)
	lon_2, v_vel_850	= PeriodicBoundaries3D(time, lon, lat, v_vel_850)
	lon_2, u_vel_250	= PeriodicBoundaries3D(time, lon, lat, u_vel_250)
	lon, v_vel_250		= PeriodicBoundaries3D(time, lon, lat, v_vel_250)

	#Determine the vorticity at 850 and 250 hPa
	vor_850			= ma.masked_all((len(time), len(lat) - 2, len(lon) - 2))
	vor_250			= ma.masked_all((len(time), len(lat) - 2, len(lon) - 2))

	for time_i in range(len(time)):
		#Vorticity for the two level
		vor_850[time_i]	= RelativeVorticity(u_vel_850[time_i], v_vel_850[time_i], x_diff, y_diff)
		vor_250[time_i]	= RelativeVorticity(u_vel_250[time_i], v_vel_250[time_i], x_diff, y_diff)
	
	#Get the correct lon/lat dimensions
	lon, lat		= lon[1:-1], lat[1:-1]

	return time, lon, lat, lat_weight, pres, U_10, temp, vor_850, vor_250

def PeriodicBoundaries2D(lon, lat, field, lon_grids = 1):
	"""Add periodic zonal boundaries for 2D field"""

	#Empty field with additional zonal boundaries
	lon_2			= np.zeros(len(lon) + lon_grids * 2)
	field_2			= ma.masked_all((len(lat), len(lon_2)))
	
	#Get the left boundary, which is the right boundary of the original field
	lon_2[:lon_grids]	= lon[-lon_grids:] - 360.0
	field_2[:, :lon_grids]	= field[:, -lon_grids:]

	#Same for the right boundary
	lon_2[-lon_grids:]	= lon[:lon_grids] + 360.0
	field_2[:, -lon_grids:]	= field[:, :lon_grids]

	#And the complete field
	lon_2[lon_grids:-lon_grids]		= lon
	field_2[:, lon_grids:-lon_grids] 	= field

	return lon_2, field_2	

def PeriodicBoundaries3D(time, lon, lat, field, lon_grids = 1):
	"""Add periodic zonal boundaries for 3D field"""

	#Empty field with additional zonal boundaries
	lon_2				= np.zeros(len(lon) + lon_grids * 2)
	field_2				= ma.masked_all((len(time), len(lat), len(lon_2)))
	
	#Get the left boundary, which is the right boundary of the original field
	lon_2[:lon_grids]		= lon[-lon_grids:] - 360.0
	field_2[:, :, :lon_grids]	= field[:, :, -lon_grids:]

	#Same for the right boundary
	lon_2[-lon_grids:]		= lon[:lon_grids] + 360.0
	field_2[:, :, -lon_grids:]	= field[:, :, :lon_grids]

	#And the complete field
	lon_2[lon_grids:-lon_grids]		= lon
	field_2[:, :, lon_grids:-lon_grids] 	= field

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
	lon_1, lat_1, lon_2, lat_2 = map(radians, [lon_1, lat_1, lon_2, lat_2]) 

	#Haversine formula 
	d_lon 	= lon_2 - lon_1 
	d_lat 	= lat_2 - lat_1 
	a 	= math.sin(d_lat/2.0)**2 + math.cos(lat_1) * math.cos(lat_2) * math.sin(d_lon/2.0)**2
	c 	= 2.0 * math.asin(sqrt(a)) 
	r 	= 6371000.0 # Radius of earth in meters
	
	return c * r #Distance between two points in meter

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

files = glob.glob(directory_data+'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'.cam2.h1.*.nc')
files.sort() #Sort the files on date 

files	= files[6:372]
print('First file: '+files[0])
print('Last file: '+files[-1])
print('#files: '+str(len(files)))

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Atmosphere_0_25_DX_DY_AREA.nc', 'r')

#Writing data to correct variable
lon_grid	= fh.variables['lon'][:]	#Longitude
lat_grid	= fh.variables['lat'][126:642]	#Latitude	
grid_x		= fh.variables['DX'][126:642]	#Zonal length of grid cell	
grid_y		= fh.variables['DY'][126:642] 	#Meridional length of grid cell		

fh.close()

lon_2, grid_x	= PeriodicBoundaries2D(lon_grid, lat_grid, grid_x)
lon, grid_y	= PeriodicBoundaries2D(lon_grid, lat_grid, grid_y)

#Generate for each grid cell the difference length
x_diff		= 0.5 * grid_x[:, 2:] + 0.5 * grid_x[:, :-2] + grid_x[:, 1:-1]
y_diff		= 0.5 * grid_y[2:] + 0.5 * grid_y[:-2] + grid_y[1:-1]

#Last elements can not be determined
x_diff	= x_diff[1:-1]
y_diff	= y_diff[:, 1:-1]

#-----------------------------------------------------------------------------------------

#Set correct ensemble member
if run_year_start == 2092:
	#Raise run ensemble by 5 (for the model output)
	run_ensemble += 5

#directory 		= '/home/rvwesten/Tropical_Cyclones/Data/RCP_'+str(run_ensemble).zfill(3)+'/'
directory 		= '/home/jasperdj/files_rene/RCP_'+str(run_ensemble).zfill(3)+'/'

if not os.path.exists(directory):
	os.mkdir(directory[:-1])
	print("Created "+directory)
if not os.path.exists(directory+'RV_Max'):
	os.mkdir(directory+'RV_Max')
	print("Created "+directory+'RV_Max') 

#-----------------------------------------------------------------------------------------

time, lon, lat, lat_weight, pres, U_10, temp, vor_850, vor_250  = ReadinData(files[0], x_diff, y_diff)

#Counters
hour_counter		= 0
RV_max_max_day		= 0

#For each day (8x3 hour) determine the pressure lows
time_RV_max		= ma.masked_all(8)
RV_max_coor_lon		= ma.masked_all((8, 800))
RV_max_coor_lat		= ma.masked_all((8, 800))
temp_anom_all		= ma.masked_all((8, 800))
#-----------------------------------------------------------------------------------------

for file_i in range(0, len(files)):
	print("File {}/{} [{}]".format(file_i+1,len(files),datetime.datetime.now()))
	print('\nProcessing '+files[file_i])
	time, lon, lat, lat_weight, pres, U_10, temp, vor_850, vor_250  = ReadinData(files[file_i], x_diff, y_diff)
		
	for time_i in range(len(time)):
		date_str = datetime.date.fromordinal(int(time[time_i]))
		#print(file_i, time_i, datetime.date.fromordinal(int(time[time_i])))
		print("{}: Time step {}/{}  [{}]".format(date_str, time_i+1, len(time),datetime.datetime.now()))

		#if os.path.exists(directory+'Atmosphere/RV_Max/RV_Max_Coordinates_'+str(datetime.date.fromordinal(int(time[time_i])))+'.nc'):
		if os.path.exists(directory+'RV_Max/RV_Max_Coordinates_'+str(datetime.date.fromordinal(int(time[time_i])))+'.nc'):
			hour_counter	= 0
			continue

		#Loop over each 3-hourly field
		lat_index_EQ	= np.where(lat >= 0.0)[0][0]

		#Get all the indices where the (absolute) vorticity is larger than 6 * 10^-5 
		index_NH	= np.where(vor_850[time_i, lat_index_EQ:] >= 6.0 * 10**(-5.0))
		lat_index_NH	= index_NH[0] + lat_index_EQ
		lon_index_NH	= index_NH[1]
		index_SH	= np.where(-vor_850[time_i, :lat_index_EQ] >= 6.0 * 10**(-5.0))
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
			vor_250_TC	= vor_250[time_i, lat_i, lon_i] * np.sign(lat[lat_i])

			if vor_850_TC < 6.0 * 10**(-5.0):
				#Relative vorticity is too low
				continue

			if (vor_850_TC - vor_250_TC) < 6.0 * 10**(-5.0):
				#No evidene for warm core
				continue

			#Get the horizontal velocity (4x4) around the eye
			lon_east	= lon_i + 14
			lon_west	= lon_i - 13
			lat_north	= lat_i + 18
			lat_south	= lat_i - 17

			if lat_south < 0:
				lat_south = 0

			#Get the dimensions for the lat field
			lat_field	= lat[lat_south:lat_north]
			weight		= lat_weight[lat_south:lat_north]
			weight		= weight / np.sum(weight)
			lon_field	= np.zeros(lon_east - lon_west)
			temp_field	= ma.masked_all((len(lat_field), len(lon_field)))
			U_10_field	= ma.masked_all((len(lat_field), len(lon_field)))

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
			lon_index_2	= (fabs(lon_field - lon[lon_i])).argmin()
			lat_index_2	= (fabs(lat_field - lat[lat_i])).argmin()
			temp_anom_TC	= temp_field[lat_index_2, lon_index_2] - temp_mean

			if temp_anom_TC < 0:
				#No warm core
				continue

			#Check whether the 10-m wind speed exceeds threshold 
			U_10_treshold = False

			for lat_j in range(len(lat_field)):
				for lon_j in range(len(lon_field)):
					#Check whether the 10-m wind speeds exceeds 10 m/s
					if U_10_field[lat_j, lon_j] < 10:
						continue

					#Check the 10-m wind speed in a 100-km search radius
					distance	= Distance(lon[lon_i], lat[lat_i], lon_field[lon_j], lat_field[lat_j]) / 1000.0

					if distance <= 100:
						#10-m wind speed exceeds treshold within 100 km
						U_10_treshold = True
						break

				if U_10_treshold:
					#Threshold is reached, stop search
					break

			if U_10_treshold == False:
				#The 10-m wind speed treshold is not valid
				continue

			#Save the coordinates and the anomalies
			RV_max_index_lon.append(lon_i)
			RV_max_index_lat.append(lat_i)
			temp_anom.append(temp_anom_TC)

		#-----------------------------------------------------------------------------------------
		distance_RV	= True
		RV_index_start	= 0

		while distance_RV:
			#Check RV maxima in the neighbourhood
			for RV_i in range(RV_index_start, len(RV_max_index_lon)):
				for RV_j in range(RV_i + 1, len(RV_max_index_lon)):
					#Get the coordinates
					lon_RV_1	= lon[RV_max_index_lon[RV_i]]
					lat_RV_1	= lat[RV_max_index_lat[RV_i]]
					lon_RV_2	= lon[RV_max_index_lon[RV_j]]
					lat_RV_2	= lat[RV_max_index_lat[RV_j]]

					if fabs(lat_RV_1 - lat_RV_2) > 3:
						#Points are too far apart
						continue

					#Determine the distance (km) between the two points
					distance	= Distance(lon_RV_1, lat_RV_1, lon_RV_2, lat_RV_2) / 1000.0

					if distance < 250:
						#Distance between pressure minima must be at least 250 km
						#Select the one with the lowest sea-level pressure
						pres_1		= pres[time_i, RV_max_index_lat[RV_i], RV_max_index_lon[RV_i]]
						pres_2		= pres[time_i, RV_max_index_lat[RV_j], RV_max_index_lon[RV_j]]

						#Retain the one with the highest RV
						index_min	= np.argmin([pres_1, pres_2])

						if index_min == 0:
							#First point has lowest PSL, remove second
							RV_max_index_lon 	= np.delete(RV_max_index_lon, RV_j)
							RV_max_index_lat 	= np.delete(RV_max_index_lat, RV_j)
							temp_anom		= np.delete(temp_anom, RV_j)

						else:
							#Second point has lowest PSL, remove first
							RV_max_index_lon 	= np.delete(RV_max_index_lon, RV_i)
							RV_max_index_lat 	= np.delete(RV_max_index_lat, RV_i)
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
					distance_RV	= False
	

		#-----------------------------------------------------------------------------------------
		#Save the indices (note that the index has a base starting at 60S)
		time_RV_max[hour_counter]				= time[time_i]
		RV_max_coor_lon[hour_counter, :len(RV_max_index_lon)]	= RV_max_index_lon
		RV_max_coor_lat[hour_counter, :len(RV_max_index_lat)]	= RV_max_index_lat
		temp_anom_all[hour_counter, :len(temp_anom)]		= temp_anom

		if len(RV_max_index_lon) > RV_max_max_day:
			#New hourly maximum of pressure minima
			RV_max_max_day = len(RV_max_index_lon)

		#Update hour counter
		hour_counter += 1

		if hour_counter == 8:
			#Remove masked arrays
			RV_max_coor_lon		= RV_max_coor_lon[:, :RV_max_max_day + 1]
			RV_max_coor_lat		= RV_max_coor_lat[:, :RV_max_max_day + 1]
			temp_anom_all		= temp_anom_all[:, :RV_max_max_day + 1]

			#Do not save the years 2002/2092 and 2008/2098
			year	= datetime.date.fromordinal(int(time[time_i])).year

			if year == 2008 or year == 2098:
				#End program
				sys.exit()

			if year > 2002 and year != 2092:
				#Save the coordinates for each day (8x 3 hour)
				#HEAT_data = netcdf.Dataset(directory+'Atmosphere/RV_Max/RV_Max_Coordinates_'+str(datetime.date.fromordinal(int(time[time_i])))+'.nc', 'w')
				HEAT_data = netcdf.Dataset(directory+'RV_Max/RV_Max_Coordinates_'+str(datetime.date.fromordinal(int(time[time_i])))+'.nc', 'w')

				HEAT_data.createDimension('time', len(time_RV_max))
				HEAT_data.createDimension('number_lows', len(RV_max_coor_lon[0]))

				HEAT_data.createVariable('time', float, ('time'), zlib=True)
				HEAT_data.createVariable('number_lows', float, ('number_lows'), zlib=True)
				HEAT_data.createVariable('lon_index', float, ('time', 'number_lows'), zlib=True)
				HEAT_data.createVariable('lat_index', float, ('time', 'number_lows'), zlib=True)
				HEAT_data.createVariable('TEMP', float, ('time', 'number_lows'), zlib=True)

				HEAT_data.variables['lon_index'].longname 	= 'Longitude index of low'
				HEAT_data.variables['lat_index'].longname 	= 'Latitude index of low (0 = 60.2S)'
				HEAT_data.variables['TEMP'].longname 		= 'Temperature anomaly w.r.t. 8x8 mean'

				HEAT_data.variables['time'].units 		= 'Days since 0001-01-01 00:00:00 UTC'
				HEAT_data.variables['TEMP'].units 		= 'deg C'

				#Writing data to correct variable	
				HEAT_data.variables['time'][:]     	= time_RV_max
				HEAT_data.variables['number_lows'][:] 	= np.arange(len(RV_max_coor_lon[0])) + 1
				HEAT_data.variables['lon_index'][:] 	= RV_max_coor_lon
				HEAT_data.variables['lat_index'][:] 	= RV_max_coor_lat
				HEAT_data.variables['TEMP'][:] 		= temp_anom_all
				
				print("Created {}".format(HEAT_data.filepath()))
				HEAT_data.close()
				

			#Reset the counters
			hour_counter	= 0
			RV_max_max_day	= 0
			
			#Make empty array for the following day
			time_RV_max		= ma.masked_all(8)
			RV_max_coor_lon		= ma.masked_all((8, 800))
			RV_max_coor_lat		= ma.masked_all((8, 800))
			temp_anom_all		= ma.masked_all((8, 800))
	
		#-----------------------------------------------------------------------------------------

