import numpy as np
import os
import utide
import matplotlib.dates as mdates
from physoce import tseries as ts
from ADCP import rditext_to_dataset
import sys

sys.path.append('..')
import datapath

# Script to read Station M ADCP text data, do tidal analysis, and export to NetCDF format

netcdf_out = datapath.adcpnc()

t1 = np.datetime64('2017-11-10T10:00Z')
t2 = np.datetime64('2018-10-17T10:00Z')    

data_path = datapath.adcptext()
data_path = os.path.join(data_path)

ds = rditext_to_dataset(data_path)

# Fill gaps 1 hour or less with linear interpolation
ds['Eas'] = ds['Eas'].interpolate_na(dim='time',limit=12)
ds['Nor'] = ds['Nor'].interpolate_na(dim='time',limit=12)

ds = ds.sel(time=slice(t1, t2))

ds['Eas_filt'] = (['time','bin'],ts.pl66(ds['Eas'],T=33*12))
ds['Nor_filt'] = (['time','bin'],ts.pl66(ds['Nor'],T=33*12))

ds.to_netcdf(netcdf_out,mode='w')