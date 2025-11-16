#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH -c 16
#SBATCH -p rome

source ~/.bashrc
conda activate geo

ulimit -v 5000000  # limit memory per process to 5GB

python -u precipitation.py ref 1 > log.ref.1 2>&1 &
python -u precipitation.py ref 2 > log.ref.2 2>&1 &
python -u precipitation.py ref 3 > log.ref.3 2>&1 &
python -u precipitation.py ref 4 > log.ref.4 2>&1 &
python -u precipitation.py ref 5 > log.ref.5 2>&1 &
python -u precipitation.py ref 6 > log.ref.6 2>&1 &

python -u precipitation.py cnt 1 > log.cnt.1 2>&1 &
python -u precipitation.py cnt 2 > log.cnt.2 2>&1 &
python -u precipitation.py cnt 3 > log.cnt.3 2>&1 &
python -u precipitation.py cnt 4 > log.cnt.4 2>&1 &
python -u precipitation.py cnt 5 > log.cnt.5 2>&1 &
python -u precipitation.py cnt 6 > log.cnt.6 2>&1 &

python -u precipitation.py sai 1 > log.sai.1 2>&1 &
python -u precipitation.py sai 2 > log.sai.2 2>&1 &
python -u precipitation.py sai 3 > log.sai.3 2>&1 &
python -u precipitation.py sai 4 > log.sai.4 2>&1 &
python -u precipitation.py sai 5 > log.sai.5 2>&1 &
python -u precipitation.py sai 6 > log.sai.6 2>&1 &
wait
