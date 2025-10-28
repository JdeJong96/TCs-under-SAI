#!/bin/bash

#SBATCH -t 48:00:00
#SBATCH -c 24
#SBATCH -p genoa
#SBATCH -o log.run-%j

source ~/.bashrc
conda activate geo

python -u Tracking_TC_RV.py ref 1 > log.ref.1 2>&1 &
python -u Tracking_TC_RV.py ref 2 > log.ref.2 2>&1 &
python -u Tracking_TC_RV.py ref 3 > log.ref.3 2>&1 &
python -u Tracking_TC_RV.py ref 4 > log.ref.4 2>&1 &
python -u Tracking_TC_RV.py ref 5 > log.ref.5 2>&1 &
python -u Tracking_TC_RV.py ref 6 > log.ref.6 2>&1 &

python -u Tracking_TC_RV.py rcp 1 > log.rcp.1 2>&1 &
python -u Tracking_TC_RV.py rcp 2 > log.rcp.2 2>&1 &
python -u Tracking_TC_RV.py rcp 3 > log.rcp.3 2>&1 &
python -u Tracking_TC_RV.py rcp 4 > log.rcp.4 2>&1 &
python -u Tracking_TC_RV.py rcp 5 > log.rcp.5 2>&1 &
python -u Tracking_TC_RV.py rcp 6 > log.rcp.6 2>&1 &

python -u Tracking_TC_RV.py sai 1 > log.sai.1 2>&1 &
python -u Tracking_TC_RV.py sai 2 > log.sai.2 2>&1 &
python -u Tracking_TC_RV.py sai 3 > log.sai.3 2>&1 &
python -u Tracking_TC_RV.py sai 4 > log.sai.4 2>&1 &
python -u Tracking_TC_RV.py sai 5 > log.sai.5 2>&1 &
python -u Tracking_TC_RV.py sai 6 > log.sai.6 2>&1 &

wait

# dump netCDF output in .ncdump files
for file in *.nc
do 
    ncdump $file > ${file}dump
done

# job statistics
seff $SLURM_JOB_ID
