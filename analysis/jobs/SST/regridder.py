""" Regrid 0.1deg POP data (tripolar grid) to the 0.25 deg CAM grid"""

import glob
import os
import numpy as np
import xarray as xr
import xesmf as xe
from load_SAIdata import Cases

# parameters for generating remap weights
# note: SRC (6th member) has a 'nice' lon field (no -1 inland), preventing regridding issues
#SRC_GRID_FILE = './data/SST.hres.sai.6.nc' # 6th member has a 'nice' lon field (no -1 inland)
SRC_GRID_FILE = ('/projects/0/nwo2021025/archive/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.006/'
                 +'ocn/hist/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.006.pop.h.2092-12.nc')
DST_GRID_FILE = ('/projects/0/nwo2021025/archive/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.001/'
                 + 'atm/hist/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.001.cam2.h0.2092-12.nc')
REMAP_METHOD  = 'conservative' # passed to xe.Regridder
WEIGHTS_FILE  = f'./data/remapweights_tx0.1v2_to_f02_t12_{REMAP_METHOD}.nc'

# append tag to output filename
DST_TAG = 'f02_t12'


def _compute_corners(ULAT, ULONG):
    """Compute grid corners. (copied from pop-tools/grid.py)"""

    nlat, nlon = ULAT.shape
    corner_lat = np.empty((nlat, nlon, 4), dtype=float)
    corner_lon = np.empty((nlat, nlon, 4), dtype=float)

    # NE corner
    corner_lat[:, :, 0] = ULAT
    corner_lon[:, :, 0] = ULONG

    # NW corner (copy from NE corner of column to the left, assume zonal periodic bc)
    corner_lat[:, :, 1] = np.roll(corner_lat[:, :, 0], 1, axis=1)
    corner_lon[:, :, 1] = np.roll(corner_lon[:, :, 0], 1, axis=1)

    # SW corner (copy from NW corner of row below, bottom row is extrapolated from 2 rows above)
    corner_lat[1:nlat, :, 2] = corner_lat[0 : nlat - 1, :, 1]
    corner_lon[1:nlat, :, 2] = corner_lon[0 : nlat - 1, :, 1]
    corner_lat[0, :, 2] = corner_lat[1, :, 2] - (corner_lat[2, :, 2] - corner_lat[1, :, 2])
    corner_lon[0, :, 2] = corner_lon[1, :, 2] - (corner_lon[2, :, 2] - corner_lon[1, :, 2])

    # SE corner (copy from NE corner of row below, bottom row is extrapolated from 2 rows above)
    corner_lat[1:nlat, :, 3] = corner_lat[0 : nlat - 1, :, 0]
    corner_lon[1:nlat, :, 3] = corner_lon[0 : nlat - 1, :, 0]
    corner_lat[0, :, 3] = corner_lat[1, :, 3] - (corner_lat[2, :, 3] - corner_lat[1, :, 3])
    corner_lon[0, :, 3] = corner_lon[1, :, 3] - (corner_lon[2, :, 3] - corner_lon[1, :, 3])

    return corner_lat, corner_lon


def gen_corner_calc(ds, cell_corner_lat_name='ULAT', cell_corner_lon_name='ULONG'):
    """
    Generates corner information and creates single dataset with output
    (adapted from https://ncar.github.io/esds/posts/2021/regrid-observations-pop-grid/)
    """

    cell_corner_lat = ds[cell_corner_lat_name]
    cell_corner_lon = ds[cell_corner_lon_name]
    ds = ds.drop_vars([cell_corner_lat_name, cell_corner_lon_name])
    
    # Use the function in pop-tools to get the grid corner information
    corn_lat, corn_lon = _compute_corners(cell_corner_lat, cell_corner_lon)

    # Make sure this returns four corner points
    assert corn_lon.shape[-1] == 4

    lon_shape, lat_shape = corn_lon[:, :, 0].shape
    out_shape = (lon_shape + 1, lat_shape + 1)

    # Generate numpy arrays to store destination lats/lons
    out_lons = np.zeros(out_shape)
    out_lats = np.zeros(out_shape)

    # Assign the northeast corner information
    out_lons[1:, 1:] = corn_lon[:, :, 0]
    out_lats[1:, 1:] = corn_lat[:, :, 0]

    # Assign the northwest corner information
    out_lons[1:, :-1] = corn_lon[:, :, 1]
    out_lats[1:, :-1] = corn_lat[:, :, 1]

    # Assign the southwest corner information
    out_lons[:-1, :-1] = corn_lon[:, :, 2]
    out_lats[:-1, :-1] = corn_lat[:, :, 2]

    # Assign the southeast corner information
    out_lons[:-1, 1:] = corn_lon[:, :, 3]
    out_lats[:-1, 1:] = corn_lat[:, :, 3]

    ds['lon_b'] = (('nlat_b', 'nlon_b'), out_lons)
    ds['lat_b'] = (('nlat_b', 'nlon_b'), out_lats)   

    return ds


def conform_dataset(ds):
    """add cell corners and rename to lat/lon for xesmf regridder"""
    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset()
    ds = gen_corner_calc(ds)
    ds = ds.rename({'TLAT':'lat', 'TLONG':'lon'})
    return ds


def gen_remap_weights():
    """generate remap weights"""
    dst = xr.open_dataset(DST_GRID_FILE)[['lat','lon']] # CAM file
    src = xr.open_dataset(SRC_GRID_FILE) # POP file
    src = conform_dataset(src)
    regrid = xe.Regridder(src, dst, method=REMAP_METHOD)
    regrid.to_netcdf(WEIGHTS_FILE)
    return


def regrid(ds, verbose=True, **kwargs):
    """regrid an xarray dataset""" 
    ds = conform_dataset(ds)
    if not os.path.exists(WEIGHTS_FILE):
        print(f'generating weights file {WEIGHTS_FILE}')
        gen_remap_weights()
    dst_grid = xr.open_dataset(DST_GRID_FILE)[['lat','lon']] # CAM file
    regrid = xe.Regridder(ds, dst_grid, method=REMAP_METHOD, reuse_weights=True,
                          weights=WEIGHTS_FILE, **kwargs)
    regridded = regrid(ds)
    return regridded


def main():
    files = sorted(glob.glob('./data/SST.hres.???.?.nc'))
    for file in files:
        ds = xr.open_dataset(file)
        regridded = regrid(ds, extrap_method='inverse_dist')
        outfile = file.rstrip('.nc') + f'.{DST_TAG}.nc'
        print(f'writing {outfile}')
        regridded.to_netcdf(outfile)
    return


if __name__ == '__main__':
    main()
