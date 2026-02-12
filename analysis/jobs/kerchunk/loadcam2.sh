#!/bin/bash

#SBATCH -t 03:00:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -o loadcam-%j

source ~/.bashrc
conda activate geo

python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.2').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.2').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.4').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.4').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.5').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.5').select('atm','h2').open_mfdataset(); ds.close()" & 
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.6').select('atm','h3').open_mfdataset(); ds.close()" &
wait
python -c "from load_SAIdata import Cases; ds=Cases('lres.sai').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('lres.sai').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('lres.sai').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('lres.sai').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('lres.cnt').select('atm','h0').open_mfdataset(); ds.close()" &
wait

