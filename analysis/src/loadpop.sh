#!/bin/bash

#SBATCH -t 03:00:00
#SBATCH -p rome
#SBATCH -c 16
#SBATCH -o loadpop-%j

source ~/.bashrc
conda activate geo

#python -c "from load_SAIdata import Cases; ds=Cases('lres.sai').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
#python -c "from load_SAIdata import Cases; ds=Cases('mres.cnt').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
#python -c "from load_SAIdata import Cases; ds=Cases('mres.sai').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &

python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.1').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.2').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.3').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.4').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.5').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.ref.6').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &

python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.1').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.2').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.3').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.4').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.5').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.cnt.6').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &

python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.1').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.2').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.3').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.4').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.5').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &
python -c "from load_SAIdata import Cases; ds=Cases('hres.sai.6').select('ocn','h.nday1').open_mfdataset_virtualizarr(ncstore_dir='~/kerchunk/pop/'); ds.close()" &

wait
