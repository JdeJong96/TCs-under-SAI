#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -o loadcam-%j

source ~/.bashrc
conda activate geo


python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.4').select('atm','h1').open_mfdataset(); ds.close()" &
wait

