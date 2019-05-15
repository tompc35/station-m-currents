import numpy as np
import os
import utide
import matplotlib.dates as mdates
from physoce import tseries as ts

# Script to read Station M ADCP text data, do tidal analysis, and export to NetCDF format

netcdf_out = 'data/MBARI_StationM_ADCP_201711_201811.nc'

t1 = np.datetime64('2017-11-10T10:00Z')
t2 = np.datetime64('2018-10-17T10:00Z')    

data_dir = '/Users/tomconnolly/work/Data/MBARI_Station_M/ADCP'
data_filename = 'MBARI_StationM_ADCP_201711_201811.txt'
data_path = os.path.join(data_dir,data_filename)

ds = rditext_to_dataset(data_path)

# Fill gaps 1 hour or less with linear interpolation
ds['Eas'] = ds['Eas'].interpolate_na(dim='time',limit=12)
ds['Nor'] = ds['Nor'].interpolate_na(dim='time',limit=12)

Eas_tide = np.nan*np.zeros(np.shape(ds['Eas']))
Nor_tide = np.nan*np.zeros(np.shape(ds['Nor']))

# tidal analysis
time = mdates.date2num(ds['time'])
for zi in np.arange(len(ds['binheight'])):
    print('bin '+str(zi+1)+'/'+str(len(ds['binheight'])))
    coef = utide.solve(time, ds['Eas'][:,zi], ds['Nor'][:,zi],
             lat=34+50/60,
             constit=['M2','S2','N2','K2','K1','O1','P1','Q1'])
    tide = utide.reconstruct(time, coef)
    Eas_tide[:,zi] = tide['u']
    Nor_tide[:,zi] = tide['v']


ds['Eas_tide'] = (['time','bin'],Eas_tide)
ds['Nor_tide'] = (['time','bin'],Nor_tide)    

ds['Eas_filt'] = (['time','bin'],ts.pl66(ds['Eas'],T=33*12))
ds['Nor_filt'] = (['time','bin'],ts.pl66(ds['Nor'],T=33*12))

ds.to_netcdf(netcdf_out)