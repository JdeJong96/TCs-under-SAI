#!/bin/bash

source ~/.bashrc
conda activate geo

# merge ensemble members
experiments=(ref cnt sai)
for exp in ${experiments[@]};
do
    cdo setmisstonn -ensmean data/SST.hres.${exp}.{1..6}.f02_t12.nc data/SST.hres.${exp}.EM1-6.nc
done

# merge experiments
python << EOF
import xarray as xr
dsets = [xr.open_dataset(f'data/SST.hres.{exp}.EM1-6.nc').assign_coords(exp=exp)
    for exp in ['ref','cnt','sai']]
xr.concat(dsets, dim='exp').to_netcdf('data/SST.nc')
EOF