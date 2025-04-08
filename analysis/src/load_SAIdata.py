import json
import os
import re
from kerchunk.netCDF3 import NetCDF3ToZarr
from kerchunk.combine import MultiZarrToZarr
import numpy as np
import dask.bag
import xarray as xr


def open_mfdataset(filepaths: list[str], ncstore_dir: str='~/kerchunk', verbose=True, **kwargs):
    """a faster alternative to xr.open_mfdataset using kerchunk
    
    This function uses kerchunk to create an NC_STORE reference file,
    that instructs the program how to read the netCDF files efficiently. 
    Coordinates must be consistent throughout all files.
    The NC_STORE is saved after first use, and will be read on each 
    subsequent usage of this function.
    
    Parameters:
    filepaths : str or list[str]
        (list of) netCDF file names, may contain wild cards
    ncstore_dir: Pathlike
        Path where NC_STORE reference files will be saved
    verbose: Bool
        Whether to print NC_STORE reference file names when reading/writing
    kwargs: dict
        any additional keyword arguments are passed on to xr.open_dataset
        
    Returns: xr.Dataset
        a Dataset instance containing all the netCDF data
        
    v0.0
    """
    
    # make sorted list of absolute filepaths
    if isinstance(filepaths, str):
        filepaths = glob.glob(filepaths)
    filepaths = sorted([os.path.abspath(fp) for fp in filepaths])
    if len(filepaths) == 1: # use xr.open_dataset directly if there is one file
        return xr.open_dataset(filepaths[0], **kwargs)
    
    # set default keyword arguments for xr.open_dataset on NC_STORE file
    required_kw = {'engine':'kerchunk', 'storage_options':{'target_protocol':'file'}}
    for (k,v) in required_kw.items():
        if k in kwargs:
            print(f'open_mfdataset(): ignoring keyword {k}')
    kwargs = kwargs | required_kw
    
    # create NC_STORE filename from netCDF filename, including timestamp
    # of first and last file. Open and return dataset if the file already exists
    ncstore_dir = os.path.expanduser(ncstore_dir)
    fparts = os.path.basename(filepaths[0]).split('.')
    fpartsf = os.path.basename(filepaths[-1]).split('.')
    fparts[-2] = f'{fparts[-2]}_{fpartsf[-2]}' # time string
    fparts[-1] = 'json' # extension
    ncstorefile = '.'.join(fparts)
    ncstore_path = os.path.join(ncstore_dir, ncstorefile)
    if not os.path.exists(ncstore_dir):
        os.mkdir(ncstore_dir)
    elif os.path.exists(ncstore_path):
        if verbose:
            print(f"Reading combined kerchunk reference file {ncstore_path}")
        return xr.open_dataset(ncstore_path, **kwargs)
    
    # make new NC_STORE data
    with xr.open_dataset(filepaths[0]) as ds:
        const_vars = [v for v in ds.variables if 'time' not in ds[v].dims]
    filebag = dask.bag.from_sequence(filepaths, npartitions=None)
    reffiles = (filebag.map(NetCDF3ToZarr, inline_threshold=0, max_chunk_size=0)
                .map(lambda z: z.translate()).compute())
    mzz = MultiZarrToZarr(reffiles, concat_dims=['time'], identical_dims=const_vars) #coo_map={'time':'cf:time'})
    
    # write NC_STORE data and return opened dataset
    with open(f"{ncstore_path}", "wb") as f:
        if verbose:
            print(f"Writing combined kerchunk reference file {ncstore_path}")
        f.write(json.dumps(mzz.translate()).encode())
    
    return xr.open_dataset(ncstore_path, **kwargs)


class Cases:
    '''Finding and opening netCDF files in all CESM1.0.4 SAI and control experiments.

    Class methods:
        select(comp, stream): select a specific model component and file stream
        open_mfdataset: open files with load_SAIdata.open_mfdataset (run after select()!)
    
    Class data:
        cases: mapping from case tags to absolute case directories
        comps: mapping from CESM component names to their module names
        files: before select(): dict of all netCDF files in case directory
               after select(): list of all netCDF files in selection
    
    Examples:
    Easy opening of datasets:
        >> from load_SAIdata import Case
        >> ds = Case(tag).select(comp,stream).open_mfdataset(*args, **kwargs)

    If tag, comp and stream are undetermined:
        >> from load_SAIdata import Case
        >> print(Case.cases)     # pick a tag
        >> mycase = Case(<tag>)
        >> print(mycase)         # pick a model component and filestream
        >> mydata = mycase.select(<comp>,<stream>)
        >> print(mydata)         # check time info from first file
        >> print(mydata.files)   # check files before opening 
        >> mydata.open_mfdataset() #
    '''

    # general info and case overview
    comps = {'atm':'cam2','ocn':'pop','lnd':'clm2','ice':'cice', 'strataero':'strataero', 'volcaero':'volcaero'} # model components
    DIR1 = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4' # root directory SAI
    DIR2 = '/projects/0/nwo2021025/archive' # root directory control
    cases = {
        'lres.spinup':f'{DIR1}/f09_g16/spinup_pd_maxcores_f09_g16/OUTPUT', # 200-300 => 2000-2100
        'lres.sai20':f'{DIR2}/lres_b.e10.B2000_CAM5.f09_g16.feedforward.001',
        'lres.sai':f'{DIR2}/lres_b.e10.B2000_CAM5.f09_g16.feedforward_2050.001',
        'mres.cnt':f'{DIR1}/f05_t12/rcp8.5_co2_f05_t12', # 2000-2100
        'mres.sai':f'{DIR2}/mres_b.e10.B2000_CAM5.f05_t12.001', # 2045-2100
        'hres.ref.1':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2002-12.001/OUTPUT',
        'hres.ref.2':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2002-12.002/OUTPUT',
        'hres.ref.3':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2002-12.003/OUTPUT',
        'hres.ref.4':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2002-12.004/OUTPUT',
        'hres.ref.5':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2002-12.005/OUTPUT',
        'hres.ref.6':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2002-12_without_SAI.001', # additional run with 6hrly 3D output
        'hres.cnt.1':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2092-12.001/OUTPUT',
        'hres.cnt.2':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2092-12.002/OUTPUT',
        'hres.cnt.3':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2092-12.003/OUTPUT',
        'hres.cnt.4':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2092-12.004/OUTPUT',
        'hres.cnt.5':f'{DIR1}/f02_t12/b.e10.B_RCP8.5_CO2_CAM5.f02_t12.started_2092-12.005/OUTPUT',
        'hres.cnt.6':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12_without_SAI.001', # additional run with 6hrly 3D output
        'hres.sai.1':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.001',
        'hres.sai.2':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.002',
        'hres.sai.3':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.003',
        'hres.sai.4':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.004',
        'hres.sai.5':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.005',
        'hres.sai.6':f'{DIR2}/hres_b.e10.B2000_CAM5.f02_t12.started_2092-12.006', # additional run with 6hrly 3D output
    }

    
    def __init__(self, tag):
        '''Initialize a specific case, identified by its tag.'''
        self.tag = tag          # tag, e.g. hres.sai.1
        self.directory = self.cases[tag] # casedirectory
        if not os.path.isdir(self.directory):
            print(f"directory does not exist: {self.directory}")
        self.name = os.path.basename(self.directory.rstrip('/OUTPUT'))  # casename
        self.model_component = None
        self.file_stream = None
        self.files = self._group_ncFiles()


    def __repr__(self):
        '''Always prints case tag and name. Additionally:
          before Case.select(): available file streams for each model component
          after Case.select(): additional time info (read first netCDF)
        '''
        msg = f'{self.tag} -- {self.name}'
        if isinstance(self.files, list): # after Case.select()
            msg += f'\nmodel component: {self.model_component}\nfile stream: {self.file_stream}'
            msg += f'\nfirst file: {self.files[0]}\nlast file: {self.files[-1]}\nnumber of files: {len(self.files)}'
            msg += self._nc_info()
        elif isinstance(self.files, dict): # before Case.select()
            msg += f'\ndirectory: {self.directory}'
            for cmp in self.files:
                msg += f'\n{cmp}: {list(self.files[cmp])}'
        return msg


    def select(self, comp, stream):
        '''Select model component and output stream'''
        self.model_component = comp
        self.file_stream = stream
        self.files = self.files[comp][stream]
        return self


    def open_mfdataset(self, *args, **kwargs):
        '''Open netCDF files, wrapper for load_SAIdata.open_mfdataset'''
        assert isinstance(self.files, list), 'attempted to open dataset without selecting a model component and file stream'
        return open_mfdataset(self.files, *args, **kwargs)


    def _nc_info(self):
        '''Open a netCDF file and return some basic info'''
        try:
            with xr.open_dataset(self.files[0]) as ds:
                msg = f'\nsteps in first file: {len(ds.time)}'
                if hasattr(ds.time, 'bounds'):
                    tbounds = ds[ds.time.bounds]
                    bnd_dim = tbounds.dims[tbounds.shape.index(2)]
                    step = tbounds.isel(time=0).diff(bnd_dim).dt.total_seconds().item()
                    msg += f'\ntime step in first file: ' + (f'{step/3600:.1f}H' if step<86400 else f'{step/86400:.1f}D' + f' ({ds.time.bounds})')
                msg += f'\ntime in first file: {[str(t) for t in ds.time.data[:3]]} ...'
        except Exception as e:
            msg = f'Could not fetch additional data from first file due to...\n{e}'
        return msg


    def _group_ncFiles(self):
        '''Find all netCDF files and group by model component and file stream.'''
        result = {comp:{} for comp in self.comps.keys()}
        for comp in self.comps:
            ncFiles = [os.path.join(root,file) 
                for root,dirs,files in os.walk(top=os.path.join(self.directory,comp)) 
                for file in files if file.endswith('.nc')]
            for file in ncFiles:
                stream = self._process_filename(file)
                if stream is None:
                    continue
                if (stream not in result[comp]) and (stream[0] in ['f','h']): # filter restarts/initial files
                    result[comp][stream] = [file]
                elif (stream in result[comp]):
                    result[comp][stream].append(file)
        for comp in result:
            result[comp] = dict(sorted(result[comp].items())) # sort streams
            for stream in result[comp]:
                result[comp][stream] = np.sort(result[comp][stream]).tolist() # sort files
        
        return result


    def _process_filename(self, fname):
        '''Read filename "fname" and return model component and file stream.'''
        if ('volcaero' in fname) or ('strataero' in fname):
            if '/copy/' in fname:
                return None, None
            stream = fname.removesuffix('.nc').split('_')[-1].removeprefix('CAM')
            if re.match('^feedback-[0-9]{4}$', stream): # group yearly files
                stream = 'feedback-YYYY'
            return stream
        parts = os.path.basename(fname).split('.')
        mods = list(self.comps.values()) # model components [cam2, pop, clm2, cice]
        try:
            imod = [(m in parts) for m in mods].index(True)
        except ValueError:
            return None, None
        is0 = parts.index(mods[imod]) + 1 # stream starts right after model component
        for is1 in range(is0,len(parts)): # some filenames have multiple stream parts, e.g. (pop).h.nday1
            if parts[is1][0].isnumeric():
                break # loop until first numeric part (date) or last part (nc)
        if re.match('^avg[0-9]{4}$', parts[is1-1]):
            parts[is1-1] = 'avgYYYY' # group yearly averages
        stream = '.'.join(parts[is0:is1])
        return stream


    def _test(self):
        '''compare found files with full directory content (using tree)'''
        print(self)
        print()
        if isinstance(self.files, dict): # before Case.select()
            for cmp in self.files:
                for strm in self.files[cmp]:
                    files = self.files[cmp][strm]
                    print(f'{".".join([cmp,strm])}: {len(files)} files')
                    print(files[0].removeprefix(self.directory))
                    if len(files)>1:
                        print(f'...\n{files[-1].removeprefix(self.directory)}')
                    print()
            os.system(f'tree {self.directory}')
        elif isinstance(self.files, list): # after Case.select()
            print(self.files[0].removeprefix(self.directory))
            if len(self.files)>1:
                print(f'...\n{self.files[-1].removeprefix(self.directory)}')
            os.system(f'tree {self.directory}/{self.model_component}')