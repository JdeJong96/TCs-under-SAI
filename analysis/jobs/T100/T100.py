# Final code for TCs-under-SAI/analysis/jobs/T100/T100.py

import os
import datetime
import flox 
import xarray as xr
from load_SAIdata import Cases
import physics
import interpolate

tslices = {  # analysis periods, inclusive
    'ref': slice('2003','2007'), 
    'cnt': slice('2093','2097'), 
    'sai': slice('2093','2097')
}

def array_job_id_to_case(i):
    """generate case tag for SLURM array job id"""
    return [case for case in Cases.cases if 'hres.' in case][int(i)]


def center_timestamps(ds):
    """center timestamps in their interval given by time_bnds"""
    # only compute the centered time from time_bnds on the first step,
    # as computing it on all steps is quite slow
    tbounds = ds[ds.time.bounds]
    if 'nbnd' in tbounds.dims:
        bnd_dim = 'nbnd'
    elif 'd2' in tbounds.dims:
        bnd_dim = 'd2'
    else:
        raise NotImplementedError(f'Bound dimension unknown. {tbounds.dims=}')
    citime = tbounds[0].mean(bnd_dim).compute() # centered initial time
    dsec = (citime - ds.time[0]).dt.total_seconds().item()
    dtshift = datetime.timedelta(seconds=dsec)
    ds['time'] = ds.time + dtshift
    return ds


def get_timestep(ds):
    """Get timestep in days"""
    dt = (ds.time[1] - ds.time[0]).dt.total_seconds().item()
    return dt/86400


def open_ds_mon_atm(exp, ens, center_time=True, **kwargs):
    """open monthly CAM dataset for one ensemble member"""
    ds = (
        Cases(f'hres.{exp}.{ens}')
        .select('atm','h0')
        .open_mfdataset(**kwargs)
        .assign_coords(exp=exp, ens=int(ens))
    )

    if center_time:
        ds = center_timestamps(ds)

    return ds


def main():
    arrjobid = os.environ['SLURM_ARRAY_TASK_ID']
    tag = array_job_id_to_case(arrjobid)
    print(f'{tag=}')
    os.makedirs('./data', exist_ok=True)
    outname = f'./data/T100.{tag}.nc'

    # read
    _,exp,ens = tag.split('.')
    ds = open_ds_mon_atm(exp, ens, chunks={})
    ds = ds.sel(time=tslices[exp])

    # interpolate temperature to 100 hPa
    P = interpolate.pressure_from_hybrid(ds)
    Ti = interpolate.interpolate(ds.T, P,  [10000], 'lev')
    Ti = Ti.rename('T100')

    # take interannual mean 
    Timean = Ti.groupby('time.month').mean()
    out = Timean.to_dataset()

    # write
    print(f'writing to {outname}')
    out.to_netcdf(outname)
    return


if __name__ == '__main__':
    main()
