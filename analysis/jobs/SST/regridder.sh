#!/bin/bash

#SBATCH -c 16
#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -o regridder-%j.out
#SBATCH -e regridder-%j.out

source ~/.bashrc
conda activate geo

python -u regridder.py

mv regridder-${SLURM_JOB_ID}.out log/
