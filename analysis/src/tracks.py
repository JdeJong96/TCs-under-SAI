import ast
import os
import glob
import numpy as np
import datetime
import xarray as xr
import nc_time_axis
import itertools
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
import shapely.geometry as sgeom
import cartopy.crs as ccrs
import cartopy.util as cutil
import cftime


def load_tracks(datapath, open_kwargs=None):
    if open_kwargs is None:
        open_kwargs = {
            'Reference': ('RCP',2002,'ref'), 
            'RCP8.5': ('RCP',2092,'rcp'), 
            'SAI2050': ('SAI',2092,'sai')
        }
    
    ds = {exp: [] for exp in open_kwargs}
    for exp, (name, year, tag) in open_kwargs.items():
        maxtn = 0
        num_days = 0
        for n in range(1,7):
            fname = f'TC_tracks.{tag}.00{n}{ext}.nc'
            fname = os.path.join(datapath, fname)
            if os.path.exists(fname):
                dsi = xr.open_dataset(fname, decode_cf=False)
                dsi = (dsi.assign_coords(id=('id',dsi.id.data+maxtn,dsi.id.attrs))
                       .assign_coords(ens=('id',np.ones(dsi.id.size)*n,{'long_name':'ensemble member'})))
                ds[exp].append(dsi)
                maxtn = dsi.id.max().item()
                times = dsi.TC_tracks.isel(data=0)
                num_days += dsi.num_days
                print(f'{tag.upper()}.00{n}: {dsi.num_days/365:.3f} years, {len(dsi.id)} tracks')
            else:
                print(f'no file named {fname}')
        ds[exp] = xr.concat(ds[exp], data_vars='minimal', dim='id', join='outer')
        ds[exp]['num_days'] = num_days.assign_attrs({'long_name':'total number of analysed days'})
    
    for k,v in ds.items():
        for i,(desc) in enumerate(v.data.description):
            shortname = desc[:desc.index(':')]
            v[shortname] = v.TC_tracks.isel(data=i).assign_attrs(
                ast.literal_eval(desc[desc.index(':')+1:]))
            v[shortname].data[v[shortname]>1e30] = np.nan
        time = v.time.copy()
        time -= 0.0625 # set to center of CESM time bounds (only for year determination)
        tmask = np.isnan(time)
        time.data = cftime.num2date(time.fillna(0), time.units, time.calendar)
        v['year'] = time.dt.year.where(~tmask, np.nan)
        v['PRECT'].data = v['PRECT']*3.6e6 # m/s to mm/hour
        v['PRECT'].attrs.update({'units':'mm/hour'})

    return ds


domains = {k:sgeom.Polygon(v) for k,v in {
    'NA': ((-65,0),(-100,20),(-100,60),(0,60),(0,0)),
    'SA': ((-65,-60),(-65,0),(25,0),(25,-60)),
    'ENP':((-180,0),(-180,60),(-100,60),(-100,20),(-65,0)),
    'WNP':((100,0),(100,60),(180,60),(180,0)),
    # 'ESP':((-180,-60),(-180,0),(-65,0),(-65,-60)),
    # 'WSP':((130,-60),(130,0),(180,0),(180,-60)),
    'SP': ((130,-60),(130,0),(295,0),(295,-60)), # last point only needed for plotting purposes
    'NI': ((40,0),(40,30),(100,30),(100,0)),
    'SI': ((25,-60),(25,0),(130,0),(130,-60)),
}.items()}


def track_stat(ds, method='min'):
    """Reduce dataset along dtime dimension through various methods.
    
    method:
        min: minimum
        max: maximum
        minmax: minimum for pressure, RMW and shear, else maximum
        mean: mean
        maxRV: maximum relative vorticity
        maxdeep: maximum 12h deepening
        develop: until maximum relative vorticity
    """
    match method:
        case 'min': # minimum along track
            out = ds.min('dtime', keep_attrs=True)
        case 'max': # maximum along track
            out = ds.max('dtime', keep_attrs=True)
        case 'minmax': # minimum for pressure/RMW/shear, else maximum
            minvars = ['PSL','PSLmon','RMW','Vshear']
            dsc = ds.copy()
            dsc[minvars] = ds[minvars].min('dtime', keep_attrs=True)
            out = dsc.max('dtime', keep_attrs=True)
        case 'mean':
            out = ds.mean('dtime', keep_attrs=True)
        case 'maxRV': # maximum relative vorticity
            out = ds.isel(dtime=ds.RV.argmax('dtime'))
        # case 'maxdeep': # max 12 hour (= 4x3hr) deepening
        #     pres = ds.pres.transpose(...,'track_time')
        #     dpres = xr.ones_like(pres[...,2:-2])*(pres.data[...,4:]-pres.data[...,:-4])
        #     times = dpres.track_time.isel(track_time=dpres.argmin('track_time'))
        #     return ds.sel(track_time=times)
        case 'maxdeep': # max 12 hour (= 4x3hr) deepening
            pres = ds.PSL.transpose(...,'dtime')
            dpres = xr.ones_like(pres[...,2:-2])*(pres.data[...,4:]-pres.data[...,:-4])
            tt = dpres.argmin('dtime')
            ttmin = (tt-2).clip(min=0)
            ttmax = (tt+2).clip(max=len(ds.dtime))
            TTmin = ds.dtime.isel(dtime=ttmin)
            TTmax = ds.dtime.isel(dtime=ttmax)
            numdays = ds.num_days
            ds = (ds.where(ds.dtime>=TTmin, np.nan)
                  .where(ds.dtime<=TTmax, np.nan)
                  .mean('dtime', keep_attrs=True))
            ds['num_days'] = numdays
            out = ds
        case 'develop': # time mean from start to max. rel. vorticity
            times = ds.dtime.isel(dtime=ds.RV.argmax('dtime'))
            numdays = ds.num_days
            ds = ds.where(ds.dtime<=times, np.nan).mean('dtime', keep_attrs=True)
            ds['num_days'] = numdays
            out = ds
    if 'dtime' in out.coords:
        out = out.reset_coords('dtime')
    out = out.set_coords('year')
    return out


def cum_prob_density(data, reverse=False):
    """Cumulative probability density of flattened data"""
    data = data.stack(x=[...]).dropna('x')
    if reverse:
        sorted = np.sort(data, axis=None)[::-1]
    else:
        sorted = np.sort(data, axis=None)
    return xr.DataArray(
        data = np.arange(data.size) / (data.size-1),
        coords = [(data.name, sorted, data.attrs)])


def prob_density(data, n=50):
    xchunks = np.array_split(np.sort(data), len(data)/n)
    for x,xs in enumerate(xchunks[:-1]):
        xchunks[x] = np.append(xchunks[x], xchunks[x+1][0])
    xfirst = [xs[0] for xs in xchunks]
    xwidths = [xs[-1]-xs[0] for xs in xchunks]
    xdiff = [(len(x)-1)/(x[-1]-x[0])/ds[exp].num_days*365 for x in xchunks]
    return xfirst, xdiff, xwidths


def prob_density(data, n=50):
    data = np.sort(data)
    x = data[n//2:-n//2]
    return x, n/(data[n:]-data[:-n])
    

# def track_density(ds, xbins=range(0,361,4), ybins=range(-90,91,4)):
#     dsi = ds[['lon','lat']]
#     if 'dtime' in dsi.dims:
#         dsi = dsi.isel(id=0).dropna('dtime', how='all')
#         hist, hlon, hlat = np.histogram2d(dsi.lon, dsi.lat, (xbins, ybins))
#         hist = hist.clip(max=1)
#         for tn in ds.id[1:]:
#             dsi = ds.sel(id=tn).dropna('dtime', how='all')
#             hist += np.histogram2d(dsi.lon, dsi.lat, (xbins, ybins))[0].clip(max=1)
#     else:
#         hist, hlon, hlat = np.histogram2d(dsi.lon, dsi.lat, (xbins, ybins))
#     hist[hist==0] = np.nan
#     hist = hist/ds.num_days.data*365
#     return hlon, hlat, hist.T


def track_density(ds, xbins=range(0,361,4), ybins=range(-90,91,4)):
    ds = ds[['lon','lat','year']]
    ds0 = ds.isel(dtime=0, drop=True, missing_dims='ignore')
    ensembles = np.unique(ds0.ens)
    years = np.unique(ds0.year)
    years = years[~np.isnan(years)]
    if (len(ensembles)>1) or (len(years)>1):
        data = []
        for ens in ensembles:
            for year in years:
                dsi = ds.where((ds0.ens==ens) & (ds0.year==year), drop=True)
                dens = track_density(dsi, xbins, ybins)
                data.append(dens.stack(x=('ens','year')))
        return xr.concat(data, 'x').unstack('x')

    if 'dtime' in ds.dims:
        dsi = ds.isel(id=0).dropna('dtime', how='all')
        hist,xbins,ybins = np.histogram2d(dsi.lon, dsi.lat, (xbins, ybins))
        hist = hist.clip(max=1)
        for tn in ds.id[1:]:
            dsi = ds.sel(id=tn).dropna('dtime', how='all')
            hist += np.histogram2d(dsi.lon, dsi.lat, (xbins, ybins))[0].clip(max=1)
    else:
        hist,xbins,ybins = np.histogram2d(ds.lon, ds.lat, (xbins, ybins))
    
    hist = xr.DataArray(
        hist,
        name = 'hcount',
        coords = {
            'xbins':('xbins',(xbins[1:]+xbins[:-1])/2,{'units':ds.lon.units}),
            'ybins':('ybins',(ybins[1:]+ybins[:-1])/2,{'units':ds.lat.units}),
        },
        dims=('xbins','ybins'),
        attrs={'long_name':'2D histogram of annual tracks'},
    ).transpose('ybins','xbins')
    hist = hist.expand_dims(('ens','year'))
    hist = hist.assign_coords({
        'ens':('ens',ds.ens.data.flatten()[0:1], ds.ens.attrs),
        'year':('year',ds.year.data.flatten()[0:1], ds.year.attrs),
    })
    xbin_edges = xr.DataArray(xbins, name='xbin_edges', 
        dims='xbin_edges', attrs={'units':ds.lon.units})
    ybin_edges = xr.DataArray(ybins, name='ybin_edges', 
        dims='ybin_edges', attrs={'units':ds.lat.units})
    return xr.Dataset({
        hist.name: hist, 
        xbin_edges.name: xbin_edges, 
        ybin_edges.name: ybin_edges
    })


def TChoursperday(ds):
    """Calculate total TC hours per day in dataset ds"""
    out = {}
    for ens in np.unique(ds.ens):
        dsi = ds.where(ds.ens==ens, drop=True)#.dropna('id',how='all')
        times_in = dsi.time.stack(x=[...]).dropna('x')
        times_out = np.arange(times_in.min(), times_in.max()+.1, 0.125)
        out[ens] = xr.DataArray(data=np.zeros(times_out.size), 
            coords={'time':('time',times_out,dsi.time.attrs)},
            name='freq', attrs={'long_name':'TC frequency','units':'TC hours per day'})
        vals, counts = np.unique(times_in, return_counts=True)
        out[ens][np.searchsorted(out[ens].time, vals)] = counts*3
        out[ens] = xr.decode_cf(out[ens].to_dataset()).freq
    out = xr.concat(out.values(), 'ens', fill_value=0)
    out['ens'] = out.ens + 1
    out = out.coarsen(time=8).sum()
    return out


def assign_dayofyear(ds):
    """Split along time into year and dayofyear dimension"""
    assert 'time' in ds.coords, "time must be a coordinate"
    assert ds.time.ndim == 1, "time must be one-dimensional"
    if not hasattr(ds.time, 'dt'):
        ds = xr.decode_cf(ds, decode_timedelta=True)
    result = []
    for year in np.unique(ds.time.dt.year):
        dsi = ds.sel(time=str(year)).drop_vars('year',errors='ignore')
        dsi = (dsi.assign_coords(time=dsi.time.dt.dayofyear)
               .rename({'time':'dayofyear'})
               .expand_dims('year'))
        result.append(dsi)
    return xr.concat(result, dim='year', join='exact')


def cyclical_rollmean(da:xr.DataArray, dim:str, w:int):
    """Apply rolling mean along cyclical dimension

    da : DataArray to apply rolling mean to
    dim : name of dimension along which to apply rolling mean
    w : window size (int)
    """
    sid = w//2
    da_x = xr.concat((
        da.isel({dim:slice(da[dim].size-sid,None)}),
        da,
        da.isel({dim:slice(None,sid)})
    ), dim=dim)
    result = da_x.rolling({dim:w}, center=True).mean(keep_attrs=True)
    result = result.isel({dim:slice(sid,da_x[dim].size-sid)})
    return result


def rolling_rolling(ds, n, rolldim, dim):
    """Expand dimensions of dataset ds by rolling the data n times
    along rolldim within a centered window of size n. Adds dimension dim.
    """
    shifts = range(n//2-n+1,n//2+1) # centered window of size n
    rolled = [ds.roll({rolldim: shift}).expand_dims(dim) for shift in shifts]
    return xr.concat(rolled, dim=dim)


def get_domain(lon, lat):
    lon = (lon + 180)%360 - 180 # convert to -180,180 range
    p = sgeom.Point(lon, lat)
    for dom in domains:
        if domains[dom].contains(p) or p.touches(domains[dom]):
            return dom
    # retry for the 180-540 range (useful for SP)
    p = sgeom.Point(lon+360, lat)
    for dom in domains:
        if domains[dom].contains(p) or p.touches(domains[dom]):
            return dom
    return None


def fix_longitude_jumps(lons, thresh=10):
    """Remove longitude jumps when crossing the Greenwich Meridian
    
    lons (1-d): longitudes
    thresh : maximum expected absolute longitude change per step
    """
    newlons = np.array(lons, copy=True)
    if any(np.abs(np.diff(lons))>thresh):
        for l in range(1,len(newlons)):
            dlon = newlons[l] - newlons[l-1]
            inc = np.sign(int(abs(dlon)>thresh) * dlon)
            if inc != 0:
                newlons[l] = newlons[l] - inc * 360
    return newlons


def track_segments(ds, clevs, x='lon', y='lat', c='PSL', reverse=True):
    """Split all tracks in ds into segments of equal intensity

    input:
    ds : xr.Dataset
        the dataset containing TC tracks
    clevs : Iterable[Numeric]
        levels (bin edges) used to discretize 'c' and segment track data
    x : str
        name of longitude variable
    y : str
        name of latitude variable
    c : str
        name of variable that will be used to segment the track data
    reverse : bool
        by default, track segments having lower values of 'c' are
        output first, this can be reversed by setting reverse=True

    returns:
    """
    ds = ds[[x,y,c]]
    cvals = 0.5 * (clevs[1:] + clevs[:-1])
    
    segments = []
    for tid in ds.id.data:
        # remove longitude jumps when crossing the Greenwich Meridian
        track = ds.sel(id=tid).dropna('dtime')
        track[x].data = fix_longitude_jumps(track[x].data)
        
        # groupby color level consecutively along track
        clev_ids = np.searchsorted(clevs, track[c])
        comb = zip(clev_ids, track[x].data, track[y].data)
        gb = itertools.groupby(comb, key=lambda x: x[0])
        track_segs = [[k,np.vstack(list(dt))[:,1:]] for k,dt in gb]
    
        # insert midway points for smooth joins of segments
        for s,(seg1, seg2) in enumerate(itertools.pairwise(track_segs)):
            ival = np.mean([seg1[1][-1],seg2[1][0]], axis=0)
            track_segs[s][1] = np.vstack([track_segs[s][1], ival])
            track_segs[s+1][1] = np.vstack([ival,track_segs[s+1][1]])
        
        segments += track_segs
    
    segments = sorted(segments, key=lambda x:x[0], reverse=reverse)
    gb = itertools.groupby(segments, key=lambda x:x[0])
    segments = {c:[segs for _,segs in grp] for c,grp in gb}
    return segments


def histogram(da, bins=np.arange(900,1015,2.5)):
    """compute histogram of da for all available ensembles and years

    input:
    da : xr.DataArray
        data to compute histogram of
    bins : Iterable[Numeric]
        bin edges used to digitize the data

    returns: xr.Dataset
        the histogram count (hcount), bin edges and central bin values
    """
    assert 'dtime' not in da.dims, "da must first be reduced, e.g. using track_stat()"
    ensembles = np.unique(da.ens)
    years = np.unique(da.year)
    data = []
    for ens in ensembles:
        ensdata = []
        for year in years:
            da_i = da.where((da.ens==ens) & (da.year==year), drop=True)
            if da_i.id.size == 0:
                continue
            hist,bins = np.histogram(da_i, bins)
            hist = xr.DataArray(
                hist,
                name = 'hcount',
                coords = {'ens':da_i.ens.isel(id=0, drop=True), 
                          'year':da_i.year.isel(id=0, drop=True), 
                          'bins':('bins',(bins[1:]+bins[:-1])/2,{'units':da.units})}, 
                dims='bins',
                attrs={'long_name':'histogram of '+getattr(da, 'long_name', da.name)},
            )
            ensdata.append(hist)
        data.append(xr.concat(ensdata, dim='year'))
    hist = xr.concat(data, join='outer', dim='ens')
    bin_edges = xr.DataArray(bins, name='bin_edges', dims='bin_edges', attrs={'units':da.units})
    return xr.Dataset({'hcount':hist, 'bin_edges':bin_edges})


def plot_histogram(ax, hists, exps=None, yunits=''):
    """make fancy histogram (difference) plot
    
    input:
    ax : matplotlib.axes.Axes
        axis to draw
    hists : dict[str,xr.Dataset]
        mapping from experiment name to Dataset containing histogram data
        Dataset is generally the output of histogram()
    exps : List[[str]
        names of experiments to plot (which must be in hists). It is possible
        to provide a list of size one with two experiments separated by a dash (-),
        in which case the difference of these experiments is plotted.
    yunits : str
        the y-axis units of a difference plot. If 'se' (not case-sensitive), the difference
        is shown in units of the combined standard error of the two datasets.

    returns:
        axis instance on which the plot is drawn
    """
    if exps is None:
        exps = list(hists.keys())
    if '-' not in exps[0]: # regular histograms
        for exp in exps:
            hist = hists[exp]
            ydata = hist.hcount.stack(x=['ens','year']).mean('x')
            ax.stairs(ydata, hist.bin_edges, fill=False, zorder=3, **default_kwargs[exp])
            binsize = (hist.bin_edges[1]-hist.bin_edges[0]).item()
            ax.set_title(f"bin={round(binsize,4)}", loc='right')
            ax.grid()
    elif '-' in exps[0]: # histogram difference plot
        for exp in exps:
            exp2, exp1 = exp.split('-')
            hist1 = hists[exp1].stack(x=['ens','year'])
            hist2 = hists[exp2].stack(x=['ens','year'])
            xdata = hist1.bins
            ydata = (hist2.hcount.mean('x') - hist1.hcount.mean('x'))
            wdata = hist1.bin_edges.data[1:] - hist1.bin_edges.data[:-1]
            stderr1 = hist1.hcount.std('x', ddof=1) / np.sqrt(hist1.hcount.x.size)
            stderr2 = hist2.hcount.std('x', ddof=1) / np.sqrt(hist2.hcount.x.size)
            stderr = np.sqrt(stderr1**2 + stderr2**2)/np.sqrt(2)
            # mask = (hist1.hcount.sum('x')>=5) & (hist2.hcount.sum('x')>=5)
            mask = (hist1.hcount.sum('x') + hist2.hcount.sum('x')) >= 10
            if yunits.lower() == 'se':
                ydata = ydata / stderr
                ax.set_ylim([-4,4])
                ax.set_yticks(range(-4,5), minor=True)
                ax.set_yticks(range(-4,5,2), minor=False)
                ax.grid(visible=True, which='minor', axis='y')
            else:
                segments = [np.column_stack([hist1.bins, stderr*i]) for i in range(-3,4)]
                ax.add_collection(LineCollection(segments, color='k', lw=0.5, zorder=4), autolim=False)
            ax.grid()
            ydata_pos = np.where((ydata >= 0) & mask, ydata, np.nan)
            ydata_neg = np.where((ydata < 0) & mask, ydata, np.nan)
            ydata_inv = np.where(~mask, ydata, np.nan)
            kwargs = default_kwargs[exp].copy()
            kwargs.pop('color')
            dh1 = ax.bar(xdata, ydata_pos, wdata, edgecolor=None, fill=True, zorder=3, color='orange', **kwargs)
            dh2 = ax.bar(xdata, ydata_neg, wdata, edgecolor=None, fill=True, zorder=3, color='lightgreen', **kwargs)
            dh3 = ax.bar(xdata, ydata, wdata, edgecolor='k', lw=1, fill=False, zorder=3)     
    return ax


def select_region(ds, region, method='maxRV'):
    """Select tracks in a specific region

    ds : xr.Dataset
        dataset containing the track data
    region : str [GLOB,NH,SH] or any key in 'domains'
        name of region
    method : str
        reduction method applied before region checking
        passed on to track_stat()

    returns : xr.Dataset
        dataset containing only the tracks that
        occur in the selected region
    """
    if 'dtime' in ds.dims:
        ds_ts = track_stat(ds, method)
        ds  = ds.assign_coords({'year': ds_ts.year})
    else:
        ds_ts = ds
    match region:
        case 'GLOB':
            return ds
        case 'NH':
            return ds.where(ds_ts.lat>=0, drop=True)
        case 'SH':
            return ds.where(ds_ts.lat<0, drop=True)
            
    ds['domain'] = ('id',[get_domain(lon, lat) for lon,lat in zip(ds_ts.lon, ds_ts.lat)])
    return ds.where(ds.domain==region, drop=True).drop_vars('domain')


def savefig(fig, fname, **kwargs):
    """calls fig.savefig(fname, **kwargs) only if fname does not exist"""
    fname = os.path.join(figpath, fname)
    if os.path.exists(fname):
        print(f"{fname} already exists, cannot overwrite")
    else:
        print(f"saving {fname}...")
        fig.savefig(fname, **kwargs)


def bootstrap(data, func, *args, n=10000, qs=(0.05,0.95), **kwargs):
    """statistical bootstrapping for significance testing
    
    data : arraylike (1D)
        data to apply bootstrapping on
    func : function
        function that returns a scalar statistic on the data, e.g. mean, median, quantile
    n : int
        number of times to resample the data
    qs : arraylike of float [0-1]
        quantile values of the statistic that are returned by this function
        with default qs, the 90% confidence interval (0.05,0.95) is returned
    args, kwargs : (keyword) arguments
        are passed on to func

    returns : np.array
        quantile values of the statistics as calculated by func on all 
        n resampled datasets
    """
    data = np.array(data)
    assert data.ndim == 1, "data must be one-dimensional"
    stat = func(data, *args, **kwargs)
    try:
        stat = np.array(stat).item()
    except ValueError:
        print(f"stat must be a scalar, but has shape {np.array(stat).shape}")
        raise
    cdata = data - stat
    ids = np.random.randint(data.size, size=(n, data.size))
    data_bs = cdata[ids] # create n samples
    stats = np.array([func(s, *args, **kwargs) for s in data_bs])
    stats += stat
    return np.quantile(stats, qs)

