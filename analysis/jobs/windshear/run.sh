#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -a 0-17

source ~/.bashrc
conda activate geo

python -u windshear.py
