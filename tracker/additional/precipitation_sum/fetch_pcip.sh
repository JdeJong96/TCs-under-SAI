#!/bin/bash

#SBATCH -t 05:00:00
#SBATCH -c 16
#SBATCH -p rome

source ~/.bashrc
conda activate geo

#python fetch_pciptrackdata.py ref 1 &
#python fetch_pciptrackdata.py ref 2 &
#python fetch_pciptrackdata.py ref 3 &
#python fetch_pciptrackdata.py ref 4 &
#python fetch_pciptrackdata.py ref 5 &
python fetch_pciptrackdata.py ref 6 &

#python fetch_pciptrackdata.py cnt 1 &
python fetch_pciptrackdata.py cnt 2 &
#python fetch_pciptrackdata.py cnt 3 &
python fetch_pciptrackdata.py cnt 4 &
python fetch_pciptrackdata.py cnt 5 &
#python fetch_pciptrackdata.py cnt 6 &

#python fetch_pciptrackdata.py sai 1 &
#python fetch_pciptrackdata.py sai 2 &
#python fetch_pciptrackdata.py sai 3 &
#python fetch_pciptrackdata.py sai 4 &
#python fetch_pciptrackdata.py sai 5 &
#python fetch_pciptrackdata.py sai 6 &

wait
