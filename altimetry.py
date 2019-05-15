import numpy as np
import pandas as pd
import xarray as xr
import gsw

def download_jpl_sla_data(date_list,lon,lat):
    '''
    Download interim gridded satellite altimetry data from:
    https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_INTERIM_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL1609
    
    Also returns geostrophic velocity calculated from sea level gradients.
    
    Note: earlier data available at:
    https://podaac.jpl.nasa.gov/dataset/SEA_SURFACE_HEIGHT_ALT_GRIDS_L4_2SATS_5DAY_6THDEG_V_JPL1609?ids=Measurement&values=Sea%20Surface%20Topography
    '''

    file_prefix = 'https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/merged_alt/L4/cdr_grid_interim/ssh_grids_v1609_'
    file_suffix = '_i.nc'

    g = 9.8

    # pre-allocate data arrays
    sla = np.nan*np.zeros(len(date_list))
    sla_err = np.nan*np.zeros(len(date_list))
    ua = np.nan*np.zeros(len(date_list))
    va = np.nan*np.zeros(len(date_list))

    print('downloading SLA data...')
    for ti,file_date in enumerate(date_list):

        print('file '+str(ti)+'/'+str(len(date_list)))

        file_datestr = (str(file_date.year) + 
                        str(file_date.month).zfill(2) +
                        str(file_date.day).zfill(2) + 
                        str(file_date.hour).zfill(2))

        file_url = file_prefix+file_datestr+file_suffix

        #with xr.open_dataset(file_url,decode_times=False) as ds_alt:
        
        for i in range(10): # try 10 times
            try:
                ds_alt = xr.open_dataset(file_url,decode_times=False)
                if ti == 0:
                    ii = int(np.argmin(np.abs(ds_alt['Latitude']-lat)))
                    jj = int(np.argmin(np.abs(ds_alt['Longitude']-lon)))
                    lat_grid = ds_alt['Latitude'][ii]
                    lon_grid = ds_alt['Longitude'][jj]
                    dlat = ds_alt['Latitude'][ii+1]-ds_alt['Latitude'][ii-1]
                    dlon = ds_alt['Longitude'][jj+1]-ds_alt['Longitude'][jj-1]
                    dy = dlat*111.32*1000 # distance in meters
                    dx = dlon*111*1000*np.cos(lat_grid*np.pi/180)
                    f = gsw.f(lat_grid) # Coriolis parameter

                sla[ti] = np.squeeze(ds_alt['SLA'][0,jj,ii])
                sla_err[ti] = np.squeeze(ds_alt['SLA_ERR'][0,jj,ii])

                ua[ti] = -(g/f)*np.squeeze(ds_alt['SLA'][0,jj,ii+1]-ds_alt['SLA'][0,jj,ii-1])/dy
                va[ti] = (g/f)*np.squeeze(ds_alt['SLA'][0,jj+1,ii]-ds_alt['SLA'][0,jj-1,ii])/dx

                ds_alt.close()
                break
            except:
                print(str(i+1)+'/10 could not open: '+file_url)
                       
    print('download finished')
            
    ds = xr.Dataset({'sla': (('time'), sla),
                 'sla_err': (('time'), sla_err),
                 'ua': (('time'), ua),
                 'va': (('time'), va)},
                {'time': date_list, 'lon': lon_grid, 'lat': lat_grid})
    
    return ds
        
if __name__ == '__main__':
    
    netcdf_out = 'data/jpl_sla_uv_timeseries.nc'
    
    t1 = np.datetime64('2015-11-26T12:00Z')
    t2 = np.datetime64('2018-06-18T12:00Z')
    date_list = pd.to_datetime(np.arange(t1,t2,np.timedelta64(5,'D')))
    
    # MBARI Station M
    lon = -122-59.9036/60+360
    lat = 35+8.4585/60
    
    ds = download_jpl_sla_data(date_list,lon,lat)
    
    ds.to_netcdf(netcdf_out)
    