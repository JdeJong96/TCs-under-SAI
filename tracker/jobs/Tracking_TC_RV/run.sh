#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -c 2
#SBATCH -p rome
#SBATCH -o log.run-%j

source ~/.bashrc
conda activate geo

python -u Tracking_TC_RV.py ref 1 > log.reg 2>&1 &
python -u Tracking_TC_RV.dev.N0.py ref 1 > log.dev.N0 2>&1 &
python -u Tracking_TC_RV.dev.N1.py ref 1 > log.dev.N1 2>&1 &
python -u Tracking_TC_RV.dev.N2.py ref 1 > log.dev.N2 2>&1 &
python -u Tracking_TC_RV.dev.N8.py ref 1 > log.dev.N8 2>&1 &
python -u Tracking_TC_RV.dev.N-1.py ref 1 > log.dev.N-1 2>&1 &
wait

# dump netCDF output in .ncdump files
for file in *.nc
do 
    ncdump $file > ${file}dump
done

# job statistics
seff $SLURM_JOB_ID
