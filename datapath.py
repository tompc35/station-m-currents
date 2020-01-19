import os.path

# Change these paths to point to the data directories
# data_dir contains subdirectories containing Station M current meter data
data_dir = '~/work/Data/MBARI_Station_M/' 
ncep_dir = '~/work/Data/NCEP-NARR/'

def rover():
    rover_csv = 'Rover_merged/Rover_II_Current_Mag_Hourly_Avg_pad_2018.csv'
    return os.path.expanduser(os.path.join(data_dir,rover_csv))

def adcpnc():
    adcp_nc = 'ADCP_netcdf/MBARI_StationM_ADCP_201711_201811.nc'
    return os.path.expanduser(os.path.join(data_dir,adcp_nc))

def adcptext():
    adcp_text = 'ADCP/MBARI_StationM_ADCP_201711_201811.txt'
    return os.path.expanduser(os.path.join(data_dir,adcp_text))

def narr():
    return os.path.expanduser(ncep_dir)