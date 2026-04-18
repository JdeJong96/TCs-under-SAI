#!/bin/bash

#SBATCH -t 01:00:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -o loadcam4-%j

source ~/.bashrc
conda activate geo


python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.1').select('atm','h0').open_mfdataset(); ds.close()" &
wait

