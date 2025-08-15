import glob
import xarray as xr

exps = ['ref','cnt','sai'] # experiments
name = 'pcip_sum' # first part of input names, also used for output name
infiles = lambda exp: f'data/{name}.{exp}.*.*.nc'
outfile = f'data/{name}.nc'

def merge_experiment(exp):
    """merge all files (for each member and year) for one experiment"""
    files = sorted(glob.glob(infiles(exp)))
    dsets = [xr.open_dataset(f).drop_vars('time') for f in files]
    ds = xr.concat(dsets, dim='x')
    ds['PRECT_std'] = ds.PRECT.std('x', ddof=1, keep_attrs=True)
    ds['Nyears'] = len(files)
    ds = ds.mean('x', keep_attrs=True)
    ds.PRECT.attrs['long_name'] = 'mean annual sum of TC precipitation'
    ds.PRECT_std.attrs['long_name'] = 'standard deviation of annual TC precipitation sum'
    ds['Nyears'].attrs['long_name'] = 'number of years the data is based on (summed over all ensemble members)' 
    return ds


def combine_all():
    """combine all files into one dataset"""
    dsets = [merge_experiment(exp) for exp in exps]
    ds = xr.concat(dsets, dim='exp').assign_coords(exp=exps)
    return ds


def main():
    ds = combine_all()
    print(f'writing to {outfile}')
    ds.to_netcdf(outfile)


if __name__ == '__main__':
    main()
