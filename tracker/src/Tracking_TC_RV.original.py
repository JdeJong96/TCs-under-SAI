#Program tracks TCs

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
#directory      = '/home/rvwesten/Tropical_Cyclones/Data/RCP/'
directory       = '/home/jasperdj/files_rene/'
run_year_start      = 2092
run_ensemble        = 1
directory_data_atm  = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'/OUTPUT/atm/hist/'
directory_data_ocn  = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'/OUTPUT/ocn/hist/daily/'

def ReadinData(filename, x_diff, y_diff):

    fh      = netcdf.Dataset(filename, 'r')

    time        = fh.variables['time'][:]       #Time
    lon     = fh.variables['lon'][:]        #Longitude
    lat     = fh.variables['lat'][126:642]      #Latitude
    lat_weight  = fh.variables['gw'][127:641]
    pres        = fh.variables['PSL'][:, 127:641] / 100.0   #Sea level pressure (hPa)
    U_10        = fh.variables['U10'][:, 127:641]   #10-meter wind speed (m/s)
    u_vel_850       = fh.variables['U850'][:, 126:642]  #Zonal velocity at 850 hPa (m/s)
    v_vel_850       = fh.variables['V850'][:, 126:642]  #Meridional velocity at 850 hPa (m/s)
    u_vel_250       = fh.variables['U250'][:, 127:641]  #Zonal velocity at 250 hPa (m/s)
    v_vel_250       = fh.variables['V250'][:, 127:641]  #Meridional velocity at 250 hPa (m/s)

    fh.close()

    #Get the precipitation file
    filename_prec   = filename[:-22] + 'h2' + filename[-20:]
    fh      = netcdf.Dataset(filename_prec, 'r')
    prec        = fh.variables['PRECT'][:, 127:641] * 1000.0 * 86400.0  #3-hourly average precipiation (mm / day)
    fh.close()

    #Set the time to toordinal form
    date        = filename[-19:-9]  
    year        = int(date[0:4])
    month       = int(date[5:7])
    day     = int(date[8:10])
    
    #Get the difference (3 hours in day)
    time_diff   = time - time[0]
    time        = datetime.datetime(year, month, day).toordinal() + time_diff + time_diff[1] / 2.0

    #Add periodic boundaries
    lon_2, u_vel_850    = PeriodicBoundaries3D(time, lon, lat, u_vel_850)
    lon, v_vel_850      = PeriodicBoundaries3D(time, lon, lat, v_vel_850)

    #Determine the vorticity at 850 and 250 hPa
    vor_850         = ma.masked_all((len(time), len(lat) - 2, len(lon) - 2))

    for time_i in range(len(time)):
        #Vorticity for the two level
        vor_850[time_i] = RelativeVorticity(u_vel_850[time_i], v_vel_850[time_i], x_diff, y_diff)
    
    #Get the correct lon/lat dimensions
    lon, lat        = lon[1:-1], lat[1:-1]

    return time, lon, lat, lat_weight, pres, U_10, vor_850, u_vel_850[:, 1:-1, 1:-1], v_vel_850[:, 1:-1, 1:-1], u_vel_250, v_vel_250, prec

def ReadinDataPSLMonth(directory_data, month, year):
    """Read in the sea-level pressure (monthly)"""

    fh = netcdf.Dataset(directory_data+'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'.cam2.h0.'+str(year)+'-'+str(month).zfill(2)+'.nc', 'r')

    pres_month      = fh.variables['PSL'][0, 127:641] / 100.0   #Sea level pressure (hPa)

    fh.close()

    return pres_month

def ReadinDataSST(directory_data, month, year):
    """Read in the Sea surface temperature (daily averages"""

    fh = netcdf.Dataset(directory_data+'b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'.pop.h.nday1.'+str(year)+'-'+str(month).zfill(2)+'-01.nc', 'r')

    time_ocn    = fh.variables['time'][:]
    time_ocn    = datetime.datetime(year, month, 1).toordinal() + time_ocn - time_ocn[0]
    lon_ocn     = fh.variables['TLONG'][:]
    lat_ocn     = fh.variables['TLAT'][:]
    KMT     = fh.variables['KMT'][:]
    SST         = fh.variables['SST'][:]    #Sea surface temperature

    fh.close()

    for time_i in range(len(time_ocn)):
        #Mask all land elements
        SST[time_i] = ma.masked_where(KMT == -1.0, SST[time_i])

    lon_ocn = ma.masked_where(KMT == -1.0, lon_ocn)
    lat_ocn = ma.masked_where(KMT == -1.0, lat_ocn)

    return time_ocn, lon_ocn, lat_ocn, SST, np.mean(SST, axis = 0)

def ReadinData_RV_Max(directory_output, time):
    """Read in the PSL minima for tracking"""

    #HEAT_data = netcdf.Dataset(directory_output+'Atmosphere/RV_Max/RV_Max_Coordinates_'+str(datetime.date.fromordinal(int(time)))+'.nc', 'r')
    HEAT_data = netcdf.Dataset(directory_output+'RV_Max/RV_Max_Coordinates_'+str(datetime.date.fromordinal(int(time)))+'.nc', 'r')

    time_RV     = HEAT_data.variables['time'][:]
    time_index  = np.where(time == time_RV)[0][0]
    lon_index   = HEAT_data.variables['lon_index'][time_index] 
    lat_index   = HEAT_data.variables['lat_index'][time_index]  
    temp_anom   = HEAT_data.variables['TEMP'][time_index]   

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
    field_2         = ma.masked_all((len(lat), len(lon_2)))
    
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
    field_2             = ma.masked_all((len(time), len(lat), len(lon_2)))
    
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
    new_field   = ma.masked_all((len(lat_field), len(lon_field)))

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
    Returns radius of maximum wind speeds (m) and the maximum wind speeds (m/s)"""

    #Get 4x4 grid around given coordinate
    lon_field, lat_field, weight, vel_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, vel_field_all)

    if reduced_radius:
        #Radius is too large, reduce region to 2 x 2
        lon_field, lat_field, weight, vel_field = FieldPartition(lon, lat, lon_index, lat_index, 7, 9, lat_weight, vel_field_all)
    
    #Find the minimum index
    index_max   = np.where(vel_field == np.max(vel_field))

    #Get the index of the local minimum of current field 
    lon_vel_max = lon_field[index_max[1][0]]
    lat_vel_max = lat_field[index_max[0][0]]

    #Determine distance (km) between pressure minima and velocity maxima
    distance    = Distance(lon[lon_index], lat[lat_index], lon_vel_max, lat_vel_max) / 1000.0

    return np.max(vel_field), distance

def VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850, v_vel_850, u_vel_250, v_vel_250):
    """Determines the vertical wind shear near pressure minimum (eye of TC)"""

    #Get 8x8 grid, centered at given coordinate
    lon_field, lat_field, weight, u_vel_850_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, u_vel_850)
    lon_field, lat_field, weight, v_vel_850_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, v_vel_850)
    lon_field, lat_field, weight, u_vel_250_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, u_vel_250)
    lon_field, lat_field, weight, v_vel_250_field = FieldPartition(lon, lat, lon_index, lat_index, 13, 17, lat_weight, v_vel_250)

    #Get weighted field
    area    = np.zeros(shape(u_vel_850_field))

    for lat_i in range(len(lat_field)):
        area[lat_i] = weight[lat_i]
    
    #Get 3 degrees from center  
    x, y    = np.meshgrid(lon_field, lat_field)
    rad_deg = np.sqrt((x - lon[lon_index])**2.0 + (y - lat[lat_index])**2.0)
    area    = ma.masked_where(np.sqrt((x - lon[lon_index])**2.0 + (y - lat[lat_index])**2.0) > 3.0, area)

    #Normalise to total
    area    = area / np.sum(area)

    #Take the spatial average
    u_vel_850_field = np.sum(u_vel_850_field * area)
    v_vel_850_field = np.sum(v_vel_850_field * area)
    u_vel_250_field = np.sum(u_vel_250_field * area)
    v_vel_250_field = np.sum(v_vel_250_field * area)
    
    u_vel_shear = fabs(u_vel_250_field - u_vel_850_field)
    v_vel_shear = fabs(v_vel_250_field - v_vel_850_field)
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

    #Retain the corresponding lat for the ocean
    index       = np.where(fabs(lat_ocn - lat_pres) < 0.2)
    index_lat   = index[0]
    index_lon   = index[1]
    lon     = lon_ocn[index]

    #Retain the corresponding lon for the ocean
    index       = np.where(fabs(lon - lon_pres) < 0.3)
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
    lon_1, lat_1, lon_2, lat_2 = map(radians, [lon_1, lat_1, lon_2, lat_2]) 

    #Haversine formula 
    d_lon   = lon_2 - lon_1 
    d_lat   = lat_2 - lat_1 
    a   = math.sin(d_lat/2.0)**2 + math.cos(lat_1) * math.cos(lat_2) * math.sin(d_lon/2.0)**2
    c   = 2.0 * math.asin(sqrt(a)) 
    r   = 6371000.0 # Radius of earth in meters
    
    return c * r #Distance between two points in meter

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------  

files = glob.glob(directory_data_atm+'3h/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_'+str(run_year_start)+'-12.'+str(run_ensemble).zfill(3)+'.cam2.h1.*-*.nc')
files.sort() #Sort the files on date
print(len(files))
files   = files[6:373]
files   = files[:2] # remove this later
#Keep track of current month (for SSTs)
current_month   = 13

print(files[0])
print(files[-1])
print(len(files))   

#-----------------------------------------------------------------------------------------

#fh = netcdf.Dataset(directory+'Land/Atmosphere_0_25_DX_DY_AREA.nc', 'r')
fh = netcdf.Dataset(directory+'Atmosphere_0_25_DX_DY_AREA.nc', 'r')

#Writing data to correct variable
lon_grid    = fh.variables['lon'][:]    #Longitude
lat_grid    = fh.variables['lat'][126:642]  #Latitude   
grid_x      = fh.variables['DX'][126:642]   #Zonal length of grid cell  
grid_y      = fh.variables['DY'][126:642]   #Meridional length of grid cell     

fh.close()

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

#-----------------------------------------------------------------------------------------

#Set correct ensemble member
if run_year_start == 2002:
    #directory_output = '/home/rvwesten/Tropical_Cyclones/Data/RCP_'+str(run_ensemble).zfill(3)+'/'
    directory_output = '/home/jasperdj/files_rene/RCP_'+str(run_ensemble).zfill(3)+'/'

else:
    #Raise run ensemble by 5 (for the model output)
    #directory_output = '/home/rvwesten/Tropical_Cyclones/Data/RCP_'+str(run_ensemble+5).zfill(3)+'/'
    directory_output = '/home/jasperdj/files_rene/RCP_'+str(run_ensemble+5).zfill(3)+'/'

#-----------------------------------------------------------------------------------------

#for file_i in range(len(files)):
for file_i in range(24):
    #Get the fields
    time, lon, lat, lat_weight, pres, U_10, vor_850, u_vel_850, v_vel_850, u_vel_250, v_vel_250, prec  = ReadinData(files[file_i], x_diff, y_diff)
    print(files[file_i])
    
    for time_i in range(len(time)):
        #loop over each 3-hourly field
        print(file_i, time_i, datetime.date.fromordinal(int(time[time_i])))
        
        #Get the coordinates for each PSL minima
        try:
            directory_RVinput = '/home/jasperdj/files_rene/RCP.started_'+str(run_year_start)+'.'+str(run_ensemble).zfill(3)
            PSL_min_lon, PSL_min_lat, temp_anom = ReadinData_RV_Max(directory_RVinput, time[time_i])
        except:
            continue

        #Retain the daily and monthly-averaged SST
        time_day    = int(datetime.date.fromordinal(int(time[time_i])).day)
        time_month  = int(datetime.date.fromordinal(int(time[time_i])).month)
        time_year   = int(datetime.date.fromordinal(int(time[time_i])).year)

        if time_year == 2008 or time_year == 2098:
            break

        if time_month != current_month:
            #New month, get the corresponding sea-level pressure and SSTs
            pres_month                  = ReadinDataPSLMonth(directory_data_atm+'monthly/', time_month, time_year)
            time_ocn, lon_ocn, lat_ocn, SST, SST_month  = ReadinDataSST(directory_data_ocn, time_month, time_year)
            current_month                   = int(datetime.date.fromordinal(int(time[time_i])).month)

        #Get the corresponding daily SST
        time_index  = np.where(time_ocn == int(time[time_i]))[0][0]
        SST_day     = SST[time_index]
        
        #Keep track of RV maxima with a sea-level pressure minima in the neighbourhood
        remove_index    = []

        for min_i in range(len(PSL_min_lat)):
            #Get the indices for the RV maxima
            lat_index   = PSL_min_lat[min_i]
            lon_index   = PSL_min_lon[min_i]

            if lat_index == 0 or lat_index == len(lat) - 1:
                #Low at 60S or 60N
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
            lon_index   = (fabs(lon_field - lon[lon_index])).argmin()
            lat_index   = (fabs(lat_field - lat[lat_index])).argmin()

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
                remove_index.append(min_i)
                continue

        #Get the indices which are still optional
        optional_index  = [i for i in range(len(PSL_min_lon)) if i not in remove_index]
        PSL_min_lon_opt = np.delete(PSL_min_lon, optional_index)
        PSL_min_lat_opt = np.delete(PSL_min_lat, optional_index)
        temp_anom_opt   = np.delete(temp_anom, optional_index)

        #Remove the false RV maxima from array
        PSL_min_lon     = np.delete(PSL_min_lon, remove_index)
        PSL_min_lat     = np.delete(PSL_min_lat, remove_index)
        temp_anom   = np.delete(temp_anom, remove_index)

        #-----------------------------------------------------------------------------------------

        for track_i in track_ID_active:
            #Loop over each active track
            exec('track = TRACK_ID_'+str(track_i))

            if len(shape(track)) == 1:
                #Only one time stamp, continue (see next step)
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
                data_track  = ma.masked_all(13)
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
                        data_track  = ma.masked_all(13)
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

        #-----------------------------------------------------------------------------------------
        #Keep track of PSL minima already assigned to active track
        remove_index    = []

        for track_i in track_ID_active:
            #Loop over each active track
            exec('track = TRACK_ID_'+str(track_i))

            if len(shape(track)) > 1:
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
                data_track  = ma.masked_all(13)
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

        #-----------------------------------------------------------------------------------------
        #Remove the already assigned PSL minima from array
        PSL_min_lon     = np.delete(PSL_min_lon, remove_index)
        PSL_min_lat     = np.delete(PSL_min_lat, remove_index)
        temp_anom   = np.delete(temp_anom, remove_index)

        for min_i in range(len(PSL_min_lat)):
            #The remaining PSL minima will start new track

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

            if SST_TC_day is ma.masked:
                #TC forms above land, discard RV maxima
                continue

            #Determine vertical wind shear
            TC_shear    = VerticalWindShear(lon, lat, lon_index, lat_index, lat_weight, u_vel_850[time_i], v_vel_850[time_i], u_vel_250[time_i], v_vel_250[time_i])

            #Find the maximum precipitation
            prec_max    = MaximumPrecipitation(lon, lat, lon_index, lat_index, prec[time_i])

            #Data for track (time, lon (pres), lat (pres), pres (3 h), pres (month), vel_max, rad. vel_max, SST (day), SST (month), vor max (1/s), TC shear (m/s), prec max (mm / day), 850 hpa temp anom (C))
            data_track  = ma.masked_all(13)
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

        #-----------------------------------------------------------------------------------------

        for track_i in track_ID_active:
            #Check the tracks whether they already finished
    
            #Get the last time stamp for each active track
            exec('track = TRACK_ID_'+str(track_i))

            if len(shape(track)) == 1:
                #Only one time stamp
                time_end_track  = track[0]
                duration_track  = 3.0
        
            else:
                time_end_track  = track[-1, 0]
                duration_track  = len(track) * 3.0

            if time[time_i] == time_end_track:
                #Still an active track, keep active
                continue

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
                        #Remove track (by setting start latitude at 90N)
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
                
                if fabs(track[0, 2]) <= 30.0 and vel_max_counter == 8 and shear_counter == 8 and track[0, 7] >= 25.0:
                    #Check whether TC originates between 30S and 30N
                    #Check whether max wind speed is sustained at 17 m/s for 24 hours
                    #Minimum sea surface temperature of 25C
                    track_ID_list.append(track_i)

                #Remove active track
                track_ID_active.remove(track_i)
            
            else:
                #Remove active track
                track_ID_active.remove(track_i)
                exec('del TRACK_ID_'+str(track_i))

        #-----------------------------------------------------------------------------------------

for track_i in track_ID_active:
    #Check the last tracks manually

    #Get the last time stamp for each active track
    exec('track = TRACK_ID_'+str(track_i))

    if len(shape(track)) == 1:
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
        
        if fabs(track[0, 2]) <= 30.0 and vel_max_counter == 8 and shear_counter == 8 and track[0, 7] >= 25.0:
            #Check whether TC originates between 30S and 30N
            #Check whether max wind speed is sustained at 17 m/s for 24 hours
            #Minimum sea surface temperature of 25C
            track_ID_list.append(track_i)

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
track_all   = ma.masked_all((len(track_ID_list), time_max, 13))

for track_i in range(len(track_ID_list)):
    #Save data to general array
    exec('track = TRACK_ID_'+str(track_ID_list[track_i]))
    track_all[track_i, :len(track)] = track

#-----------------------------------------------------------------------------------------
#Save the coordinates for each track
#HEAT_data = netcdf.Dataset(directory_output+'Atmosphere/TC_tracker_RV.nc', 'w')
HEAT_data = netcdf.Dataset(directory_output+'TC_tracker_RV.nc', 'w')

HEAT_data.createDimension('track_number', len(track_ID_list))
HEAT_data.createDimension('time', time_max)
HEAT_data.createDimension('data', 13)

HEAT_data.createVariable('track_number', float, ('track_number'), zlib=True)
HEAT_data.createVariable('time', float, ('time'), zlib=True)
HEAT_data.createVariable('data', float, ('data'), zlib=True)
HEAT_data.createVariable('TC_tracks', float, ('track_number', 'time', 'data'), zlib=True)

HEAT_data.variables['track_number'].longname    = 'Track number'
HEAT_data.variables['time'].longname        = 'Lifetime TC'
HEAT_data.variables['data'].longname        = 'Time, lon (pres), lat (pres), pres (3-h), pres (month), max vel, rad. max vel, SST (day), SST (month), RV max (1/s), wind shear (m/s), prec max (mm / day)'

HEAT_data.variables['time'].units       = 'hours'

#Writing data to correct variable   
HEAT_data.variables['track_number'][:]  = np.arange(len(track_ID_list)) + 1
HEAT_data.variables['time'][:]      = np.arange(time_max) * 3
HEAT_data.variables['data'][:]      = np.arange(13)
HEAT_data.variables['TC_tracks'][:]     = track_all

HEAT_data.close()


