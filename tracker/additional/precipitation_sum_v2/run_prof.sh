#!/bin/bash

#SBATCH -t 04:30:00
#SBATCH -c 1
#SBATCH -p rome

source ~/.bashrc
conda activate geo

python -u precipitation.py
