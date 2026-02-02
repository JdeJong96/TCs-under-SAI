#!/bin/bash

#SBATCH -c 16
#SBATCH -t 00:10:00
#SBATCH -p rome
#SBATCH -o landmask-%j.out
#SBATCH -e landmask-%j.out

source ~/.bashrc
conda activate geo

python -u << EOF
"""Generate a land mask and regrid to the CAM grid"""
import xarray as xr
import regridder

OUTFILE = './data/landmask.f02_t12.nc'

KMT = xr.open_dataset(regridder.SRC_GRID_FILE, decode_timedelta=False).KMT 
KMT = KMT.astype('int') # num ocean levels
mask = (KMT <= 0).astype(KMT.TLONG.dtype) # 0: ocean, 1: land
gausmask = regridder.regrid(mask, ignore_degenerate=True)
gausmask.to_netcdf(OUTFILE)
EOF

rm -rf ./__pycache__
mv landmask-${SLURM_JOB_ID}.out log/