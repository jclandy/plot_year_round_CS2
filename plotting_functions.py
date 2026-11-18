"""Python counterparts of the supplied UiT sea-ice MATLAB workflows.

Time is always axis 0; spatial arrays follow the latitude variable's dimensions.
"""
from pathlib import Path
from datetime import datetime, timedelta
import re
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from scipy.io import loadmat
from scipy.spatial import cKDTree
from scipy import ndimage

MASIE = ['Beaufort Sea', 'Chukchi Sea', 'East Siberian Sea', 'Laptev Sea',
         'Kara Sea', 'Barents Sea', 'Greenland Sea', 'Baffin Bay/Gulf of St. Lawrence',
         'Canadian Archipelago', 'Hudson Bay', 'Central Arctic', 'Bering Sea',
         'Baltic Sea', 'Sea of Okhotsk', 'Yellow Sea', 'Cook Inlet']
ASSETS = Path(__file__).parent / 'aux_files'
# Change these defaults once to adjust every exported figure.
FONT_SIZE = 14
MAP_TITLE_SIZE = 17
LINE_COLORS = dict(snow='#7b3294', freeboard='#008837', thickness='#ca0020',
                   concentration='#182dff', myi='#ffc457')
VARIABLES = dict(ifb='sea_ice_freeboard', ifb_unc='sea_ice_freeboard_unc',
                 sit='sea_ice_thickness', sit_unc='sea_ice_thickness_unc',
                 sd='snow_depth', sd_unc='snow_depth_unc',
                 sic='sea_ice_concentration', type='sea_ice_type')


def _crs(a, e, phi_c, lambda_0):
    from pyproj import CRS
    return CRS.from_proj4(f'+proj=stere +lat_0={90 if phi_c > 0 else -90} '
                         f'+lat_ts={phi_c} +lon_0={lambda_0} +a={a} '
                         f'+b={a*np.sqrt(1-e*e)} +units=m')


def polarstereo_fwd(phi, lambda_, a=6378137., e=0.08181919, phi_c=-70., lambda_0=0.):
    """Latitude/longitude degrees to metres; defaults match supplied MATLAB helper."""
    from pyproj import Transformer
    return Transformer.from_crs('EPSG:4326', _crs(a,e,phi_c,lambda_0), always_xy=True).transform(lambda_, phi)


def polarstereo_inv(x, y, a=6378137., e=0.08181919, phi_c=-70., lambda_0=0.):
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(_crs(a,e,phi_c,lambda_0), 'EPSG:4326', always_xy=True).transform(x,y)
    return lat, lon


def _project(lat, lon):
    return polarstereo_fwd(lat, lon, phi_c=70.)


def viridis(n=256):
    return plt.get_cmap('viridis')(np.linspace(0,1,n))[:, :3]


def inpaint_nans(a, method=1):
    """Nearest valid neighbour fill (replacement, not MATLAB method-1 parity)."""
    if method != 1:
        raise ValueError('Only replacement method=1 is supported')
    a = np.asarray(a, float).copy()
    bad = ~np.isfinite(a)
    if bad.any() and not bad.all():
        idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
        a[bad] = a[tuple(idx[:, bad])]
    return a


def _mean(a, axis=0):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(a, axis=axis)


class Archive:
    """Read each file eagerly and close it; no dask requirement or open file accumulation."""
    def __init__(self, folder, version):
        import xarray as xr
        records = []
        for p in Path(folder).glob(f'uit_sit_*_{version}.nc'):
            match = re.search(r'_(\d{8})_', p.name)
            if match:
                records.append((datetime.strptime(match[1], '%Y%m%d'), p))
        records.sort()
        if not records:
            raise FileNotFoundError(f'No uit_sit_*_{version}.nc files in {folder}')
        self.dates, self.files = zip(*records)
        if len(set(self.dates)) != len(self.dates):
            raise ValueError('Duplicate dates: select a directory containing one product only')
        self.tvec = np.array([[d.year,d.month,d.day] for d in self.dates])
        self.slots = self.tvec[:,1]*2 + (self.tvec[:,2] > 15)
        with xr.open_dataset(self.files[0], decode_times=False) as ds:
            self.dims = ds.latitude.squeeze().dims
            self.lat = ds.latitude.squeeze().values.astype(float)
            self.lon = ds.longitude.squeeze().transpose(*self.dims).values.astype(float)
        if self.lat.ndim != 2:
            raise ValueError('Expected two-dimensional latitude/longitude grids')
        self.x,self.y = _project(self.lat,self.lon)
        self.points = np.column_stack((self.x.ravel(),self.y.ravel()))
        if not np.isfinite(self.points).all():
            raise ValueError('Grid coordinates must be finite')
        self.tree = cKDTree(self.points)

    def read(self, variables, indices=None):
        import xarray as xr
        result = {k:[] for k in variables}
        for p in self.files:
            with xr.open_dataset(p, decode_times=False) as ds:
                if not np.allclose(ds.latitude.squeeze().transpose(*self.dims),self.lat,equal_nan=True) or not np.allclose(ds.longitude.squeeze().transpose(*self.dims),self.lon,equal_nan=True):
                    raise ValueError(f'Grid changes in {p}')
                for key,var in variables.items():
                    a = ds[var]
                    for dim in list(a.dims):
                        if dim not in self.dims:
                            if a.sizes[dim] != 1:
                                raise ValueError(f'{var}: unexpected non-singleton dimension {dim}')
                            a = a.isel({dim:0})
                    values = a.transpose(*self.dims).values.astype(float)
                    result[key].append(values if indices is None else values.ravel()[indices])
        return {k:np.stack(v) for k,v in result.items()}


def _clim(a, slots, minimum=1):
    means = {}
    for slot in np.unique(slots):
        selected = a[slots==slot]
        mean = np.asarray(_mean(selected)).copy()
        mean = np.where(np.isfinite(selected).sum(axis=0)>=minimum,mean,np.nan)
        means[slot] = mean
    return means, a-np.stack([means[s] for s in slots])


def _out(folder, version):
    path=Path(folder)/f'general_plots_{version}'
    path.mkdir(parents=True,exist_ok=True)
    return path


def _safe(name):
    return re.sub(r'[^\w .-]+','_',name)


def _series(archive, idx, name, outfolder, plot_opt, version, region=False,
            geometry=None, coastlines=True):
    idx=np.unique(np.asarray(idx,dtype=int))
    if not len(idx):
        raise ValueError('Selection contains no grid centres; increase radius or change region')
    values=archive.read(VARIABLES,idx)
    values['sid']=values['sit']-values['ifb']
    missing=np.isnan(values['sid']) & (np.isnan(values['sic']) | (values['sic']<15))
    values['sid'][missing]=0
    # Deliberately preserve the source formula, not independent-error propagation.
    values['sid_unc']=values['sit_unc']-values['ifb_unc']
    data={'selection_is_region':region,'tvec':archive.tvec, 'dates':np.array(archive.dates,dtype='datetime64[ns]'),
          'lat_idx':archive.lat.ravel()[idx], 'lon_idx':archive.lon.ravel()[idx],
          'coverage_frac_t':np.isfinite(values['sit']).mean(axis=1)}
    for k,v in values.items():
        data[k+'_t_idx']=v
        data[k+'_t']=_mean(v,axis=1)
    if region:
        bad=(data['coverage_frac_t']<.05)|(data['ifb_t']<0)
        for k in values:
            if k not in ('sic','type'):
                data[k+'_t'][bad]=0
    for k in values:
        clim,anom=_clim(data[k+'_t'],archive.slots)
        data[k+'_clim']=np.array([clim[s] for s in sorted(clim)])
        if not k.endswith('_unc'):
            data[k+'_t_anom']=anom
    data['climatology_slots']=np.array(sorted(set(archive.slots)))
    if geometry:
        data.update(geometry)
    if plot_opt:
        _plot_series(data,name,_out(outfolder,version)/f'time_series_{_safe(name)}_{version}.png', coastlines=coastlines)
    return data


def _polar_inset(fig, slot, d, coastlines=True, *, rect=None):
    import cartopy.crs as ccrs
    projection = ccrs.NorthPolarStereo(central_longitude=-45)
    ax = (fig.add_axes(rect, projection=projection) if rect is not None
          else fig.add_subplot(slot, projection=projection))
    lat = np.asarray(d['lat_idx'])
    lower = max(45., min(80., float(np.nanmin(lat))-6))
    # Set a symmetric projected extent: geographic corner bounds can clip a
    # rotated polar map and exclude the selected location.
    px,py=projection.transform_point(-45.,lower,ccrs.PlateCarree())
    radius=float(np.hypot(px,py))
    ax.set_xlim(-radius,radius); ax.set_ylim(-radius,radius)
    ax.set_facecolor('white')
    if coastlines:
        import cartopy.feature as cfeature
        ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='#cecce3',zorder=1)
        ax.coastlines(resolution='50m', linewidth=.4,color='#666666',zorder=2)
    ax.scatter(d['lon_idx'], lat, s=5, color='#6baed6',
               transform=ccrs.PlateCarree(), zorder=3)
    if 'lat_gate' in d:
        ax.plot(d['lon_gate'], d['lat_gate'], color='#08519c', linewidth=2,
                transform=ccrs.PlateCarree(), zorder=4)
    if 'polygon_lat' in d:
        ax.plot(np.r_[d['polygon_lon'],d['polygon_lon'][0]],
                np.r_[d['polygon_lat'],d['polygon_lat'][0]], color='#08519c',
                linewidth=1.5, transform=ccrs.PlateCarree(), zorder=4)
    if 'reference_lat' in d:
        ax.scatter(d['reference_lon'],d['reference_lat'],s=45,color='#08519c',
                   transform=ccrs.PlateCarree(),zorder=4)
    # Zoom the square inset around the selected cells and supplied geometry.
    all_lat=list(np.asarray(d['lat_idx']).ravel())
    all_lon=list(np.asarray(d['lon_idx']).ravel())
    for la,lo in [('lat_gate','lon_gate'),('polygon_lat','polygon_lon'),('reference_lat','reference_lon')]:
        if la in d:
            all_lat.extend(np.atleast_1d(d[la])); all_lon.extend(np.atleast_1d(d[lo]))
    xy=projection.transform_points(ccrs.PlateCarree(),np.asarray(all_lon),np.asarray(all_lat))
    centre=(xy[:,:2].min(axis=0)+xy[:,:2].max(axis=0))/2
    half=max(np.ptp(xy[:,:2],axis=0).max()*.75,450000.)
    ax.set_xlim(centre[0]-half,centre[0]+half)
    ax.set_ylim(centre[1]-half,centre[1]+half)

    return ax


def _plot_series(d,name,path,*,coastlines=True):
    with plt.rc_context({'font.size': FONT_SIZE, 'axes.labelsize': FONT_SIZE,
                         'legend.fontsize': FONT_SIZE-1}):
        fig=plt.figure(figsize=(15,13))
        gs=fig.add_gridspec(4,1,left=.09,right=.98,bottom=.07,top=.96,hspace=.12)
        axs=[]
        for i in range(4):
            axs.append(fig.add_subplot(gs[i,0],sharex=axs[0] if axs else None))
        t=d['dates']; a=axs[0]
        a.fill_between(t,-d['sid_t'],d['ifb_t'],color='#83a7cb',alpha=1,linewidth=0,label='Sea ice thickness')
        a.fill_between(t,d['ifb_t'],d['ifb_t']+d['sd_t'],color='#b5d8eb',alpha=1,linewidth=0,label='Snow thickness')
        a.set(title=name.replace('_',' '),ylabel='Relative height ASL [m]',ylim=(-4,2))
        axs[1].plot(t,d['sd_t_anom'],color=LINE_COLORS['snow'],lw=2,label='Snow depth')
        axs[1].plot(t,d['ifb_t_anom'],color=LINE_COLORS['freeboard'],lw=2,label='Sea ice freeboard')
        axs[1].set_ylabel('Thickness anomaly [m]')
        axs[2].plot(t,d['sit_t_anom'],color=LINE_COLORS['thickness'],lw=2,label='Sea ice thickness')
        axs[2].set_ylabel('Thickness anomaly [m]')
        if d.get('selection_is_region',False):
            values=np.r_[d['sd_t_anom'],d['ifb_t_anom']]
            bound=round(float(np.nanmax(np.abs(values))),1)+.05 if np.isfinite(values).any() else .3
            thickness=d['sit_t_anom']
            tb=round(float(np.nanmax(np.abs(thickness))))+.5 if np.isfinite(thickness).any() else 2.
            axs[1].set_ylim(-bound,bound); axs[2].set_ylim(-tb,tb)
        else:
            axs[1].set_ylim(-.3,.3); axs[2].set_ylim(-2,2)
        axs[3].plot(t,d['sic_t'],color=LINE_COLORS['concentration'],lw=2,label='Sea ice concentration')
        axs[3].plot(t,d['type_t']*100,color=LINE_COLORS['myi'],lw=2,label='Fraction of MYI')
        axs[3].set(ylabel='Percentage [%]',ylim=(0,100),xlabel='Year')
        from matplotlib.dates import YearLocator, DateFormatter
        axs[0].set_title(name.replace('_',' '),fontweight='bold')
        for i,a in enumerate(axs):
            a.grid(alpha=.25)
            a.legend(loc='lower right' if i==3 else 'upper right',fancybox=False,framealpha=1)
            a.tick_params(labelbottom=i==3,direction='in',top=True,right=True,length=7)
            a.set_xlim(t[0],t[-1])
        axs[-1].xaxis.set_major_locator(YearLocator())
        axs[-1].xaxis.set_major_formatter(DateFormatter('%Y'))
        fig.canvas.draw()
        box=axs[-1].get_position()
        _polar_inset(fig,None,d,coastlines,rect=[box.x0+.01,box.y0+.018,.10,.10])
        fig.savefig(path,dpi=180,bbox_inches='tight',pad_inches=.08); plt.close(fig)


def time_series_point(name, lat, lon, sit_folder, search_radius, outfolder, plot_opt, version, *, coastlines=True):
    if search_radius < 0: raise ValueError('search_radius must be nonnegative metres')
    a=Archive(sit_folder,version)
    idx=a.tree.query_ball_point(_project(lat,lon),search_radius)
    return _series(a,idx,name,outfolder,plot_opt,version,
                   geometry=dict(reference_lat=lat,reference_lon=lon),coastlines=coastlines)


def time_series_region(name, lat, lon, masie_region, sit_folder, outfolder, plot_opt, version,
                       *, mask_path=ASSETS/'arctic_region_mask_MASIE.mat', coastlines=True):
    a=Archive(sit_folder,version)
    geometry=None
    if masie_region is not None and np.size(masie_region):
        region=int(masie_region)
        if region not in range(1,17): raise ValueError('MASIE region must be 1..16')
        mask=loadmat(mask_path)
        x,y=_project(mask['mask_lat'],mask['mask_lon'])
        valid=np.isfinite(x)&np.isfinite(y)
        tree=cKDTree(np.column_stack((x[valid],y[valid])))
        nearest=tree.query(a.points)[1]
        idx=np.flatnonzero(mask['region_mask'][valid][nearest]==region)
        name=MASIE[region-1]
    else:
        lat=np.asarray(lat,float); lon=np.asarray(lon,float)
        if lat.size==2 and lon.size==2:
            lat=lat[[0,0,1,1]]; lon=lon[[0,1,1,0]]
        if lat.size<3 or lat.shape!=lon.shape:
            raise ValueError('Supply matching polygon vertices or two latitude/longitude bounds')
        geometry=dict(polygon_lat=lat,polygon_lon=lon)
        x,y=_project(lat,lon)
        polygon=MplPath(np.column_stack((x,y)))
        idx=np.flatnonzero(polygon.contains_points(a.points,radius=1e-6))
    return _series(a,idx,name,outfolder,plot_opt,version,region=True,
                   geometry=geometry,coastlines=coastlines)


def _gate_points(shp_name, *, gate_crs=None):
    """Decode .prj or explicit CRS; detect geographic coordinates with stale projected .prj."""
    import shapefile
    from pyproj import CRS, Transformer, Geod
    prj=Path(shp_name).with_suffix('.prj')
    declared=CRS.from_user_input(gate_crs) if gate_crs is not None else (CRS.from_wkt(prj.read_text()) if prj.exists() else None)
    parts=[]
    geod=Geod(a=6378137,b=6378137)
    with shapefile.Reader(str(shp_name)) as reader:
        for shape in reader.shapes():
            if shape.shapeType not in (3,13,23):
                raise ValueError('Gate must contain polyline features')
            boundaries=list(shape.parts)+[len(shape.points)]
            for start,end in zip(boundaries[:-1],boundaries[1:]):
                points=np.asarray(shape.points[start:end],dtype=float)[:,:2]
                if len(points)<2 or not np.isfinite(points).all():
                    raise ValueError('Gate needs at least two finite vertices')
                geographic=np.all(np.abs(points[:,0])<=180) and np.all(np.abs(points[:,1])<=90)
                source=declared
                if gate_crs is None and declared is not None and declared.is_projected and geographic:
                    warnings.warn('Gate coordinates look like longitude/latitude but .prj declares a projected CRS. '
                                  'Using EPSG:4326; supply gate_crs explicitly to override.',UserWarning,stacklevel=2)
                    source=CRS.from_epsg(4326)
                if source is None:
                    if not geographic:
                        raise ValueError('Projected gate requires .prj or gate_crs')
                    source=CRS.from_epsg(4326)
                transform=Transformer.from_crs(source,'EPSG:4326',always_xy=True)
                lo,la=transform.transform(points[:,0],points[:,1])
                points=np.column_stack((lo,la))
                if not np.isfinite(points).all() or np.any(np.abs(points[:,1])>90):
                    raise ValueError('Invalid transformed gate coordinates; check gate_crs')
                dense=[tuple(points[0])]
                for p,q in zip(points[:-1],points[1:]):
                    dist=geod.inv(*p,*q)[2]
                    if dist<1e-6: continue
                    n=max(0,int(np.ceil(dist/1113.1949))-1)
                    if n: dense.extend(geod.npts(*p,*q,n))
                    dense.append(tuple(q))
                parts.append(np.array(dense))
    if not parts: raise ValueError('Empty gate')
    if len(parts)!=1: raise ValueError('Use one connected polyline part per gate')
    if len(parts[0])<2: raise ValueError('Gate has zero length')
    return parts[0][:,1],parts[0][:,0]


def _gate_direction(lat, lon):
    """Describe increasing distance using the gate's first and last vertices."""
    from pyproj import Geod
    azimuth, _, distance = Geod(ellps='WGS84').inv(
        float(lon[0]), float(lat[0]), float(lon[-1]), float(lat[-1]))
    if distance < 1:
        return 'Start to end (vertex order)'
    # Classify the overall endpoint bearing; wrap correctly across the dateline.
    if abs(np.sin(np.deg2rad(azimuth))) >= abs(np.cos(np.deg2rad(azimuth))):
        return 'West to East' if azimuth > 0 else 'East to West'
    return 'South to North' if abs(azimuth) < 90 else 'North to South'


def _plot_gate(d,name,path,*,coastlines=True):
    with plt.rc_context({'font.size':FONT_SIZE}):
        fig=plt.figure(figsize=(15,13))
        gs=fig.add_gridspec(4,2,width_ratios=[30,.5],left=.09,right=.93,
                            bottom=.07,top=.96,hspace=.12,wspace=.025)
        axs=[]
        labels={'ifb':'Sea ice freeboard [m]','sd':'Snow depth [m]',
                'sit':'Sea ice thickness [m]','type':'Fraction of MYI [%]'}
        for i,(k,limit) in enumerate(zip(['ifb','sd','sit','type'],[.5,.5,4,100])):
            ax=fig.add_subplot(gs[i,0],sharex=axs[0] if axs else None); axs.append(ax)
            z=d[k+'_gate_t']*(100 if k=='type' else 1)
            from matplotlib.colors import LinearSegmentedColormap
            if k=='type':
                cm=LinearSegmentedColormap.from_list('myi_brown',['#f5f5f5','#f6e8c3','#dfc27d','#bf812d','#8c510a'])
            else:
                cm=plt.get_cmap('PuBu' if k in ('ifb','sd') else 'YlGnBu').copy()
            cm.set_bad('#ffffff')
            image=ax.pcolormesh(d['dates'],d['D_gate']/1000,np.ma.masked_invalid(z),
                                shading='nearest',vmin=0,vmax=limit,cmap=cm)
            ax.set_ylabel('Distance along gate [km]')
            ax.set_ylim(d['D_gate'][0]/1000,d['D_gate'][-1]/1000)
            ax.tick_params(labelbottom=i==3)
            if not np.isfinite(z).any():
                ax.text(.5,.5,'No valid data in this gate buffer',ha='center',transform=ax.transAxes)
            fig.colorbar(image,cax=fig.add_subplot(gs[i,1]),label=labels[k])
        fig.canvas.draw()
        box=axs[-1].get_position()
        inset = _polar_inset(fig,None,d,coastlines,rect=[box.x0+.01,box.y0+.07,.10,.10])
        direction = _gate_direction(d['lat_gate'], d['lon_gate'])
        inset.text(.5, -.12, 'Distance increases\n' + direction,
                   transform=inset.transAxes, ha='center', va='top',
                   fontsize=FONT_SIZE-3, clip_on=False,
                   bbox=dict(facecolor='white',edgecolor='none',alpha=.9))
        from matplotlib.dates import YearLocator,DateFormatter
        axs[-1].xaxis.set_major_locator(YearLocator())
        axs[-1].xaxis.set_major_formatter(DateFormatter('%Y'))
        for ax in axs:
            ax.tick_params(direction='in',top=True,right=True,length=7)
        axs[0].set_title(name,fontweight='bold'); axs[-1].set_xlabel('Year')
        fig.savefig(path,dpi=180,bbox_inches='tight',pad_inches=.08); plt.close(fig)


def time_series_gate_shp(shp_name,sit_folder,search_radius,outfolder,plot_opt,version,
                         *,gate_crs=None,coastlines=True):
    if search_radius<0: raise ValueError('search_radius must be nonnegative metres')
    lat,lon=_gate_points(shp_name,gate_crs=gate_crs)
    a=Archive(sit_folder,version); x,y=_project(lat,lon)
    groups=a.tree.query_ball_point(np.column_stack((x,y)),search_radius)
    idx=np.unique([i for group in groups for i in group]).astype(int)
    name=Path(shp_name).stem
    geometry=dict(lat_gate=lat,lon_gate=lon,D_gate=np.r_[0,np.cumsum(np.hypot(np.diff(x),np.diff(y)))])
    d=_series(a,idx,name,outfolder,plot_opt,version,geometry=geometry,coastlines=coastlines)
    lookup={i:j for j,i in enumerate(idx)}
    for k in (*VARIABLES,'sid','sid_unc'):
        v=d[k+'_t_idx']
        d[k+'_gate_t']=np.stack([_mean(v[:,[lookup[i] for i in g]],axis=1) if g else np.full(len(a.files),np.nan) for g in groups])
    d['gate_cell_count']=np.array([len(g) for g in groups])
    if plot_opt:
        _plot_gate(d,name,_out(outfolder,version)/f'time_series_across_gate_{_safe(name)}_{version}.png',coastlines=coastlines)
    return d

# (range, anomaly bound, standard-deviation bound, colormap, units)
STYLES = {
 'radar_freeboard':((0,.35),.15,.1,'Spectral_r','m'),
 'sea_ice_freeboard':((0,.5),.2,.15,'Spectral_r','m'),
 'sea_ice_thickness':((0,4),1.5,1,'Spectral_r','m'),
 'snow_depth':((0,.5),.2,.15,'Blues','m'),
 'sea_ice_concentration':((0,100),30,30,'PuBu_r','%'),
 'sea_ice_type':((0,1),.3,.3,'Blues',''),
 'sea_ice_density':((868,918),24,10,'YlGnBu','kg/m³')}


def _clean(field,sic,lat,*,zero_missing=True):
    v=field.copy()
    filled=inpaint_nans(v)
    hole=(lat>88)&np.isnan(v)
    v[hole]=filled[hole]
    v[(lat>88.5)|(sic<15)]=np.nan
    if zero_missing: v[v==0]=np.nan
    labels,n=ndimage.label(np.isfinite(v),structure=np.ones((3,3)))
    sizes=np.bincount(labels.ravel())
    v[(labels>0)&(sizes[labels]<2)]=np.nan
    return v


def _map_panel(a,panels,parameter,path,*,columns=4,coastlines=True):
    """Square polar panels styled after the supplied MATLAB examples."""
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    if not panels: raise ValueError('No observations in requested period')
    limits,anomaly,std,cmap,unit=STYLES[parameter]
    rows=int(np.ceil(len(panels)/columns))
    projection=ccrs.NorthPolarStereo(central_longitude=0,true_scale_latitude=70)
    kinds=list(dict.fromkeys(panel[3] for panel in panels))
    # Symmetric square extent, north pole at centre; corners extend south of 68 N.
    _,edge=projection.transform_point(0,68,ccrs.PlateCarree())
    radius=abs(edge)
    with plt.rc_context({'font.size':FONT_SIZE}):
        monthly=len(kinds)>1 and 'anomaly' in kinds
        panel=3.15; gap_x=.16; gap_y=.20; margin=.18
        bottom=.90 if monthly else margin
        width=2*margin+columns*panel+(columns-1)*gap_x+(0 if monthly else 1.0)
        height=margin+bottom+rows*panel+(rows-1)*gap_y
        fig=plt.figure(figsize=(width,height))
        map_axes=[]
        last_ax=None
        for i,(title,field,sic,kind) in enumerate(panels):
            row,col=divmod(i,columns)
            left=(margin+col*(panel+gap_x))/width
            low=(bottom+(rows-1-row)*(panel+gap_y))/height
            ax=fig.add_axes([left,low,panel/width,panel/height],projection=projection)
            map_axes.append(ax)
            last_ax=ax
            ax.set_xlim(-radius,radius); ax.set_ylim(-radius,radius)
            ax.set_facecolor('white')
            lo,hi=limits if kind=='value' else ((-anomaly,anomaly) if kind=='anomaly' else (0,std))
            cm=cmap if kind=='value' else ('RdBu' if kind=='anomaly' else 'viridis')
            v=_clean(field,sic,a.lat,zero_missing=kind!='anomaly')
            if np.isfinite(v).any():
                ax.contourf(a.x,a.y,np.ma.masked_invalid(v),levels=np.linspace(lo,hi,65),
                            cmap=cm,extend='both',corner_mask=True,transform=projection)
            if np.isfinite(sic).any() and np.nanmin(sic)<15<np.nanmax(sic):
                ax.contour(a.x,a.y,sic,levels=[15],colors='black',linewidths=.5,transform=projection)
            if coastlines:
                # Finer coastlines and lavender land match MATLAB's map background.
                ax.add_feature(cfeature.LAND.with_scale('50m'),facecolor='#cecce3',
                               edgecolor='none',zorder=4)
                ax.coastlines(resolution='50m',linewidth=.35,color='#666666',zorder=5)
            ax.text(.995,.01,title,transform=ax.transAxes,ha='right',va='bottom',
                    fontsize=MAP_TITLE_SIZE,fontweight='bold',zorder=6,
                    bbox=dict(facecolor='white',edgecolor='none',pad=.4))
            ax.spines['geo'].set_linewidth(.4)
        fig.canvas.draw()
        # Independent vertical bars in the right margin, aligned with the last row.
        box=last_ax.get_position()
        for i,kind in enumerate(kinds):
            lo,hi=limits if kind=='value' else ((-anomaly,anomaly) if kind=='anomaly' else (0,std))
            cm=cmap if kind=='value' else ('RdBu' if kind=='anomaly' else 'viridis')
            label=('Standard deviation' if kind=='std' else parameter.replace('_',' ').capitalize()+(' anomaly' if kind=='anomaly' else ''))
            barunit='percentage points' if kind=='anomaly' and unit=='%' else unit
            if len(kinds)>1 and 'anomaly' in kinds:
                # Monthly comparisons: two separate horizontal bars below last pair.
                # Anchor each bar to the actual final-row map, with a 0.12-inch gap.
                target=map_axes[-min(columns,len(map_axes))+max(0,columns-2)+i].get_position()
                cax=fig.add_axes([target.x0,target.y0-.28/height,target.width,.16/height])
                orientation='horizontal'
            else:
                # Seasonal means and SD have their own row-aligned vertical bars.
                row_box=(map_axes[i*columns+columns-1].get_position() if len(kinds)>1 and rows==2 else box)
                cax=fig.add_axes([row_box.x1+.12/width,row_box.y0,.16/width,row_box.height])
                orientation='vertical'
            fig.colorbar(ScalarMappable(norm=Normalize(lo,hi),cmap=cm),cax=cax,
                         orientation=orientation,ticks=np.linspace(lo,hi,5),
                         label=label+(f' [{barunit}]' if barunit else ''))
            cax.tick_params(labelsize=FONT_SIZE)
            cax.xaxis.label.set_fontweight('bold'); cax.yaxis.label.set_fontweight('bold')
            if orientation=='horizontal' and kind=='anomaly':
                cax.tick_params(axis='x',labelrotation=45)
        fig.savefig(path,dpi=180,bbox_inches='tight',pad_inches=.08); plt.close(fig)
    return path


def _interval_title(archive,index):
    """Use each file's actual time centre and interval length when supplied."""
    import xarray as xr
    with xr.open_dataset(archive.files[index]) as ds:
        if 'time' not in ds or 'time_interval' not in ds:
            return archive.dates[index].strftime('%b %d')
        value=np.asarray(ds['time'].values).reshape(-1)[0]
        if np.issubdtype(np.asarray(value).dtype,np.datetime64):
            centre=value.astype('datetime64[us]').astype(datetime)
        elif isinstance(value,(int,float,np.integer,np.floating)):
            centre=datetime(1970,1,1)+timedelta(seconds=float(value))
        else:
            centre=datetime(value.year,value.month,value.day)
        interval=np.asarray(ds['time_interval'].values).reshape(-1)[0]
        if np.issubdtype(np.asarray(interval).dtype,np.timedelta64):
            days=float(interval/np.timedelta64(1,'D'))
        else:
            days=float(interval)
            units=ds['time_interval'].attrs.get('units','days').lower()
            if units.startswith('second'): days/=86400
            elif units.startswith('hour'): days/=24
        half=timedelta(days=int(np.floor(days/2)))
        return (centre-half).strftime('%b %d')+' - '+(centre+half).strftime('%b %d')


def _map_data(parameter,sit_folder,version):
    if parameter not in STYLES: raise ValueError(f'Unknown parameter: {parameter}')
    a=Archive(sit_folder,version)
    v=a.read({'p':parameter,'sic':'sea_ice_concentration'})
    return a,v['p'],v['sic']


def plot_map_climatology(parameter,sit_folder,outfolder,version,*,coastlines=True):
    """Seasonal mean and sample standard deviation, as in supplied script."""
    a,p,sic=_map_data(parameter,sit_folder,version)
    panels=[]
    for months,label in [([12,1,2],'Dec, Jan, Feb'),([3,4,5],'Mar, Apr, May'),
                         ([6,7,8],'Jun, Jul, Aug'),([9,10,11],'Sep, Oct, Nov')]:
        selected=np.isin(a.tvec[:,1],months)
        mean=_mean(p[selected]); concentration=_mean(sic[selected])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            std=np.nanstd(p[selected],axis=0,ddof=1)
        panels.extend([(label,mean,concentration,'value'),(label,std,concentration,'std')])
    path=_out(outfolder,version)/f'{parameter}_seasonal_climatology_{version}.png'
    panels=panels[::2]+panels[1::2]
    return _map_panel(a,panels,parameter,path,columns=4,coastlines=coastlines)


def plot_map_year(year,parameter,sit_folder,outfolder,version,*,minimum_count=5,coastlines=True):
    a,p,sic=_map_data(parameter,sit_folder,version)
    means,_=_clim(p,a.slots,minimum_count)
    anom=np.stack([v-means[s] for v,s in zip(p,a.slots)])
    indices=np.flatnonzero(a.tvec[:,0]==year)
    if not len(indices): raise ValueError(f'No data for {year}')
    paths=[]
    for values,kind,infix in [(p,'value',''),(anom,'anomaly','anomalies_')]:
        panels=[(_interval_title(a,i),values[i],sic[i],kind) for i in indices]
        path=_out(outfolder,version)/f'{parameter}_{infix}{year}_{version}.png'
        paths.append(_map_panel(a,panels,parameter,path,coastlines=coastlines))
    return paths


def plot_map_month(month,parameter,sit_folder,outfolder,version,*,minimum_count=5,coastlines=True):
    if month not in range(1,13): raise ValueError('month must be 1..12')
    a,p,sic=_map_data(parameter,sit_folder,version)
    means,_=_clim(p,a.slots,minimum_count)
    anomalies=np.stack([v-means[s] for v,s in zip(p,a.slots)])
    panels=[]
    for year in np.unique(a.tvec[a.tvec[:,1]==month,0]):
        select=(a.tvec[:,0]==year)&(a.tvec[:,1]==month)
        concentration=_mean(sic[select])
        panels.extend([(str(year),_mean(p[select]),concentration,'value'),
                       (str(year),_mean(anomalies[select]),concentration,'anomaly')])
    path=_out(outfolder,version)/f'{parameter}_month_{month}_{version}.png'
    return _map_panel(a,panels,parameter,path,coastlines=coastlines)


