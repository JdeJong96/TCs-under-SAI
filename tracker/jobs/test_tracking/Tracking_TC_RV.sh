#!/bin/bash
#SBATCH -c 16  
#SBATCH -p rome
#SBATCH -t 2:30:00 # takes about 50-80 minutes per model year
#SBATCH -o Tracking_TC_RV.log%j

echo "START: $(date)"
source ~/.bashrc
conda activate geo
export OMP_NUM_THREADS=1 # prevent subprocesses from using paralellism
#python -u Tracking_TC_RV.py ref 2 > log.reg.ref 2>&1 &
python -u Tracking_TC_RV.dev.py ref 1 > log.seeds.ref 2>&1 &
#python -u Tracking_TC_RV.py sai 1 > log.reg.sai 2>&1 &
#python -u Tracking_TC_RV.dev.py sai 1 > log.dev.sai 2>&1 &
wait
echo "END: $(date)"
seff $SLURM_JOB_ID

