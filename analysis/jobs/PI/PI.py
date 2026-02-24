import os
import datetime
import flox 
import numpy as np
import xarray as xr
import tcpyPI
from load_SAIdata import Cases
import physics
import interpolate

NAME = 'PI'
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
    outname = f'./data/{NAME}.{tag}.nc'

    # read
    _,exp,ens = tag.split('.')
    ds = open_ds_mon_atm(exp, ens, chunks={})
    ds = ds.sel(time=TSLICES[exp])

    # average to 0.5 degrees (more manageable results)
    ds = ds.coarsen(lat=2, lon=2).mean()

    # interpolate to pressure levels
    p = interpolate.pressure_from_hybrid(ds)/100
    plev = np.append(interpolate.PLEV/100, [1010,1020,1030,1040])
    dsp = interpolate.interpolate(ds[['T','Q','PSL','TS','OCNFRAC']], 
            p, plev, dim='lev', on_missing_core_dim='copy')
    dsp = dsp.isel(plev=slice(None,None,-1))
    
    # unit conversion
    dsp['T'] = (dsp.T-273.15).assign_attrs({'units':'°C'})
    dsp['Q'] = (dsp.Q*1000).assign_attrs({'units':'g/kg'})
    dsp['TS'] = (dsp.TS - 273.15).assign_attrs({'units':'°C'})
    dsp['PSL'] = (dsp.PSL/100).assign_attrs({'units':'hPa'})

    # compute interannual mean
    dsp = dsp.groupby('time.month').mean().compute()
    
    # compute potential intensity
    vmax, pmin, ifl, t0, otl = xr.apply_ufunc(
        tcpyPI.pi, dsp.TS, dsp.PSL, dsp.plev, dsp.T, dsp.Q,
        kwargs={'miss_handle':0},
        input_core_dims=[[],[],['plev'],['plev'],['plev']],
        output_core_dims=[[],[],[],[],[]],
        output_dtypes=[dsp.T.dtype]*5,
        vectorize=True,
    )
    
    # add metadata
    dsp['vmax'] = vmax.assign_attrs({'long_name':'Maximum Potential Intensity', 'units':'m/s'})
    #dsp['pmin'] = pmin.assign_attrs({'long_name':'Minimum Central Pressure', 'units':'hPa'})
    #dsp['ifl'] = ifl.assign_attrs({'long_name':'pyPI Flag'})
    dsp['t0'] = t0.assign_attrs({'long_name':'Outflow Temperature', 'units':'K'})
    dsp['otl'] = otl.assign_attrs({'long_name':'Outflow Temperature Level', 'units':'hPa'})

    # mask and reduce file size (only TS needed for further analysis)
    dsp = dsp.where(dsp.OCNFRAC>=0.9)
    dsp = dsp.drop_vars(['T','Q','PSL'])

    # compute logs of individual contributing terms
    lnpi, lneff, lndiseq, lnCKCD = xr.apply_ufunc(
        tcpyPI.decompose_pi,
        dsp.vmax,
        dsp.TS+273.15,
        dsp.t0,
        input_core_dims=[[],[],[]],
        output_core_dims=[[],[],[],[]],
        output_dtypes=[dsp.TS.dtype]*4,
        vectorize=True,    
    )
    
    # add metadata
    dsp['lnpi'] = lnpi.assign_attrs({'long_name':'Natural log(Potential Intensity)'})
    dsp['lneff'] = lneff.assign_attrs({'long_name':'Natural log(Tropical Cyclone Efficiency)'})
    dsp['lndiseq'] = lndiseq.assign_attrs({'long_name':'Natural log(Thermodynamic Disequilibrium)'})
    dsp['lnCKCD'] = lnCKCD.assign_attrs({'long_name':'Natural log(Ck/CD)','units':'unitless constant'})

    # add constants and weights
    dsp['lnCKCD'] = dsp.lnCKCD.mean(...)
    dsp['gw'] = ds.gw

    # write
    print(f'writing to {outname}')
    dsp.to_netcdf(outname)
    return


if __name__ == '__main__':
    main()
