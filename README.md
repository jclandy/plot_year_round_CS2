# Basic plotting and time series extraction routines for the UiT CryoSat-2 Year-round Arctic sea ice thickness dataset

## Data
[Year-round Arctic sea ice thickness from CryoSat-2 Baseline-E observations 2010-2025, V2](https://doi.org/10.21334/NPOLAR.2026.6581F94E)

## MATLAB Example

### Installation
[othercolor](https://se.mathworks.com/matlabcentral/fileexchange/30564-othercolor)

### Plotting
Point to input directory holding data and output directory for plots 
```
sit_folder = '/.../year_round_sit_v2p1'
parameter = 'sea_ice_thickness'
outfolder = cd
version = 'v2p1'
```

Create biweekly sea-ice thickness maps for 2020 and for all Septembers between 2011 and 2024
```
plot_map_year(2020,parameter,sit_folder,outfolder,version)
plot_map_month(9,parameter,sit_folder,outfolder,version)
```
Create plot of climatological sea-ice thickness and its variability between years
```
plot_map_climatology(parameter,sit_folder,outfolder,version)
```
Create plots and extract time series data for 80N 0E, for the Central Arctic MASIE region, and for Fram Strait
```
point = time_series_point('Reference', 80, 0, sit_folder, 100e3, outfolder, 1, version)
region = time_series_region('', [], [], 11, sit_folder, outfolder, 1, version)
gate = time_series_gate_shp('FRAM.shp',sit_folder,100e3,outfolder,1,version)
```

## Python Example

### Installation
```python
conda create -n year_round_cs2 -c conda-forge python=3.11 numpy scipy matplotlib xarray netcdf4 pyproj pyshp cartopy
conda activate year_round_cs2
```

### Plotting
Point to input directory holding data and output directory for plots 
```python
from pathlib import Path
import numpy as np
from plotting_functions import (plot_map_climatology, plot_map_year, plot_map_month,
                    time_series_point, time_series_region, time_series_gate_shp)
                    
folder = Path('/.../year_round_sit_v2p1')
out = Path('/.../')
version = 'v2p1'
assets = Path(__file__).parent / 'aux_files'
```

Create biweekly sea-ice thickness maps for 2020 and for all Septembers between 2011 and 2024
```python
plot_map_year(2020, 'sea_ice_thickness', folder, out, version)
plot_map_month(9, 'sea_ice_thickness', folder, out, version)
```
Create plot of climatological sea-ice thickness and its variability between years
```python
plot_map_climatology('sea_ice_thickness', folder, out, version)
```
Create plots and extract time series data for 80N 0E, for the Central Arctic MASIE region, and for Fram Strait
```python
point = time_series_point('Reference', 80., 0., folder, 100e3, out, 1, version)
region = time_series_region('', None, None, 11, folder, out, 1, version)
# np.savez_compressed(out/'timeseries_central_arctic.npz', **region)
gate = time_series_gate_shp(assets / 'FRAM.shp', folder, 100e3, out, 1, version, gate_crs='EPSG:4326')
```
