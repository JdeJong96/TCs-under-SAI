#!/bin/bash
#SBATCH -c 24
#SBATCH -p genoa
#SBATCH -t 40:00:00 # takes about 50-80 minutes per model year
#SBATCH -o Tracking_TC_RV.seeds.log%j

source ~/.bashrc
conda activate geo
export OMP_NUM_THREADS=1 # prevent subprocesses from using paralellism
script=Tracking_TC_RV.py
logname=log.seeds

#python -u $script ref 1 > $logname.ref.1 2>&1 &
#python -u $script ref 2 > $logname.ref.2 2>&1 &
#python -u $script ref 3 > $logname.ref.3 2>&1 &
#python -u $script ref 4 > $logname.ref.4 2>&1 &
#python -u $script ref 5 > $logname.ref.5 2>&1 &
python -u $script ref 6 > $logname.ref.6 2>&1 &

#python -u $script rcp 1 > $logname.rcp.1 2>&1 &
#python -u $script rcp 2 > $logname.rcp.2 2>&1 &
#python -u $script rcp 3 > $logname.rcp.3 2>&1 &
#python -u $script rcp 4 > $logname.rcp.4 2>&1 &
#python -u $script rcp 5 > $logname.rcp.5 2>&1 &
python -u $script rcp 6 > $logname.rcp.6 2>&1 &

python -u $script sai 1 > $logname.sai.1 2>&1 &
python -u $script sai 2 > $logname.sai.2 2>&1 &
python -u $script sai 3 > $logname.sai.3 2>&1 &
python -u $script sai 4 > $logname.sai.4 2>&1 &
python -u $script sai 5 > $logname.sai.5 2>&1 &
python -u $script sai 6 > $logname.sai.6 2>&1 &
wait

seff $SLURM_JOB_ID
