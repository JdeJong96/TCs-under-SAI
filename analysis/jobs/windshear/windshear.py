import os
import datetime
import flox 
import xarray as xr
from load_SAIdata import Cases
import physics


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
    arrjobid = os.environ['SLURM_ARRAY_TASK_ID']
    tag = array_job_id_to_case(arrjobid)
    print(f'{tag=}')
    _,exp,ens = tag.split('.')
    ds = open_ds_3h_atm(exp, ens)
    ds = ds.sel(time=tslices[exp])

    # monthly mean VWS
    VWS = physics.windshear_250_850(ds)
    VWSmonmean = VWS.resample(time='1MS').mean().groupby('time.month').mean()
    
    # monthly mean days with low VWS
    dt = get_timestep(ds)
    VWSdays = (VWS <= 12.5) * dt
    VWSdays = VWSdays.resample(time='1MS').sum().groupby('time.month').mean()
    VWSdays = VWSdays.astype('float32').assign_attrs({
        'long_name': 'number of days having at most 12.5 m/s wind shear'})

    outname = f'../../data/VWS.{tag}.nc'
    out = xr.Dataset({'VWS': VWSmonmean, 'lowVWSdays': VWSdays})
    print(f'writing to {outname}')
    out.to_netcdf(outname)
    return


if __name__ == '__main__':
    main()
