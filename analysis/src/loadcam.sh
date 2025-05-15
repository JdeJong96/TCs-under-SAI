#!/bin/bash

#SBATCH -t 03:00:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -o loadcam-%j

source ~/.bashrc
conda activate geo

python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.1').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.1').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.1').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.1').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.1').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.1').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.1').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.1').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.1').select('atm','h5').open_mfdataset(); ds.close()" &
wait
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.2').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.2').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.2').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.2').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.2').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.2').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.2').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.2').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.2').select('atm','h5').open_mfdataset(); ds.close()" &
wait
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.3').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.3').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.3').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.3').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.3').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.3').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.3').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.3').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.3').select('atm','h5').open_mfdataset(); ds.close()" &
wait
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.4').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.4').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.4').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.4').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.4').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.4').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.4').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.4').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.4').select('atm','h5').open_mfdataset(); ds.close()" &
wait
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.5').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.5').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.5').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.5').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.5').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.5').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.5').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.5').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.5').select('atm','h5').open_mfdataset(); ds.close()" &
wait
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.6').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.6').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.6').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.6').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.6').select('atm','h1').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.6').select('atm','h2').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.6').select('atm','h3').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.6').select('atm','h4').open_mfdataset(); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.6').select('atm','h5').open_mfdataset(); ds.close()" &
wait

