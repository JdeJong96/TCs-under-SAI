#!/bin/bash
#SBATCH -c 16  
#SBATCH -p genoa
#SBATCH -t 14:00:00 # takes about 50-80 minutes per model year
#SBATCH -o Tracking_TC_RV.seeds.log%j

source ~/.bashrc
conda activate geo
export OMP_NUM_THREADS=1 # prevent subprocesses from using paralellism
script=Tracking_TC_RV
python -u $script.py ref 1 > $script.log.ref.1 &
python -u $script.py ref 2 > $script.log.ref.2 &
python -u $script.py ref 3 > $script.log.ref.3 &
python -u $script.py ref 4 > $script.log.ref.4 &
python -u $script.py ref 5 > $script.log.ref.5 &
python -u $script.py ref 6 > $script.log.ref.6 &

#python -u $script.py rcp 1 > $script.log.rcp.1 &
#python -u $script.py rcp 2 > $script.log.rcp.2 &
#python -u $script.py rcp 3 > $script.log.rcp.3 &
#python -u $script.py rcp 4 > $script.log.rcp.4 &
#python -u $script.py rcp 5 > $script.log.rcp.5 &
python -u $script.py rcp 6 > $script.log.rcp.6 &

python -u $script.py sai 1 > $script.log.sai.1 &
python -u $script.py sai 2 > $script.log.sai.2 &
python -u $script.py sai 3 > $script.log.sai.3 &
python -u $script.py sai 4 > $script.log.sai.4 &
python -u $script.py sai 5 > $script.log.sai.5 &
python -u $script.py sai 6 > $script.log.sai.6 &
wait
