#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -a 0-17

source ~/.bashrc
conda activate geo

python -u Uzm.py

mkdir -p log
mv slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out log/
