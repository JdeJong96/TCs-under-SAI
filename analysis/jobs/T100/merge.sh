#!/bin/bash

#SBATCH -c 16
#SBATCH -p rome
#SBATCH -t 00:01:00

source ~/.bashrc
conda activate geo

python << EOF
import xarray as xr

# take ensemble mean
for exp in ['ref','cnt','sai']:
    with xr.open_mfdataset(f'data/T100.hres.{exp}.?.nc', combine='nested', concat_dim='ens') as ds:
        ds.mean('ens',keep_attrs=True).to_netcdf(f'data/T100.hres.{exp}.EM1-6.nc')

# merge experiments
dsets = [xr.open_dataset(f'data/T100.hres.{exp}.EM1-6.nc').assign_coords(exp=exp)
    for exp in ['ref','cnt','sai']]
xr.concat(dsets, dim='exp').to_netcdf('data/T100.nc')
EOF
