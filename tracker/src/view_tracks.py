import netCDF4 as netcdf
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

directory = '/home/jasperdj/files_rene/RCP_001/'
fh1 = netcdf.Dataset(directory+'RV_Max/RV_Max_Coordinates_2003-01-10.nc','r')
fh2 = netcdf.Dataset(directory+'TC_tracker_RV.nc','r')
print([var for var in fh1.variables])
print([var for var in fh2.variables])
fh2vars = ['track_number','time','data','TC_tracks']

lons = fh2.variables['TC_tracks'][:,:,1]
lats = fh2.variables['TC_tracks'][:,:,2]

fig = plt.figure(figsize=(12,8))
ax = fig.add_axes([0.125,0.125,0.875,0.875], projection = ccrs.PlateCarree())
for track_id in range(len(lons)):
	ax.plot(lons[track_id,:], lats[track_id,:])
plt.savefig('CycloneTracks.png')
plt.show()

fh1.close()
fh2.close()