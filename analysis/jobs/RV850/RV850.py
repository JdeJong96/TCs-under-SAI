import os
import datetime
import flox 
import xarray as xr
from load_SAIdata import Cases
import physics


NAME = 'RV850'
TSLICES = {  # analysis periods, inclusive
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
    citime = ds.time_bnds[0].mean('nbnd').compute() # centered initial time
    dsec = (citime - ds.time[0]).dt.total_seconds().item()
    dtshift = datetime.timedelta(seconds=dsec)
    ds['time'] = ds.time + dtshift
    return ds


def get_timestep(ds):
    """Get timestep in days"""
    dt = (ds.time[1] - ds.time[0]).dt.total_seconds().item()
    return dt/86400


def open_ds_3h_atm(exp, ens, center_time=True, **kwargs):
    """open 3-hourly atmospheric dataset for one ensemble member"""
    # file stream corresponding to 3-hourly data
    if (exp in ['ref','cnt']) and (int(ens) in range(1,6)):
        stream = 'h1'
    else:
        stream = 'h5'

    ds = (
        Cases(f'hres.{exp}.{ens}')
        .select('atm',stream)
        .open_mfdataset(**kwargs)
        .assign_coords(exp=exp, ens=int(ens))
    )

    if center_time:
        ds = center_timestamps(ds)

    return ds


def main():
    # select dataset
    arrjobid = os.environ['SLURM_ARRAY_TASK_ID']
    tag = array_job_id_to_case(arrjobid)
    print(f'{tag=}')

    # set output
    os.makedirs('./data', exist_ok=True)
    outname = f'./data/{NAME}.{tag}.nc'

    # open dataset
    _,exp,ens = tag.split('.')
    ds = open_ds_3h_atm(exp, ens)
    ds = ds.sel(time=TSLICES[exp])

    # monthly mean RV computation
    RV850 = physics.relative_vorticity_850(ds)
    RV850 = RV850.rename(NAME)
    RV850mon = RV850.resample(time='1MS').mean().groupby('time.month').mean()
    out = RV850mon.to_dataset()
    out = out.compute()

    # write
    print(f'writing to {outname}')
    out.to_netcdf(outname)

    return


if __name__ == '__main__':
    main()
