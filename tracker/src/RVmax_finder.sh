#!/bin/bash
#SBATCH -c 192 # adjust depending on # of unfinished tasks (dry run RV_max_finder)
#SBATCH -p rome,genoa
#SBATCH -t 2:00:00  # processing 1 year of data takes 52min on 64 cores
#SBATCH -o RVmax_finder.seeds.log.ref-%j

source ~/.bashrc
conda activate geo
export OMP_NUM_THREADS=1 # prevent subprocesses from using paralellism
cd /home/jasperdj/files_rene
echo "experiment: $1"
echo "run ensemble: $2"
python -u RVmax_finder.seeds.py $1 $2
