#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -c 2
#SBATCH -p rome
#SBATCH -o log.run-%j

source ~/.bashrc
conda activate geo

python -u Tracking_TC_RV.py ref 1 > log.reg 2>&1 &
python -u Tracking_TC_RV.dev.py ref 1 > log.dev 2>&1 &

wait

seff $SLURM_JOB_ID
