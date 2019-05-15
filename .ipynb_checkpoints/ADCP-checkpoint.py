import numpy as np
import pandas as pd
import xarray as xr

def rditext_to_dataset(data_path,time_span=None):
    '''
    Load RDI text data file into xarray dataset.
    
    INPUTS:
    data_path: path to text file created by RDI software.
    time_span (optional): list of two datetime64 values defining time span of data to load.
    
    RETURNS
    xarray dataset
    
    EXAMPLE:
    import numpy as np
    data_path = 'RDI_ADCP_data.txt'
    t1 = np.datetime64('2017-11-10T10:00Z')
    t2 = np.datetime64('2018-10-17T10:00Z')  
    ds = rditext_to_dataset(data_path,[t1,t2])
    '''
    
    nhead = 16
    headerlines = list()

    with open(data_path) as f: 
        for i in range(nhead):
            line = f.readline()
            line = line.rstrip('\n')
            line = line.replace('\"','')
            headerlines.append(line)

    # create list of variable names from header line
    var_names = headerlines[12].split('\t')

    df = pd.read_csv(data_path,
                     skip_blank_lines=False,
                     header=None,
                     names=var_names,
                     skiprows=nhead,
                     delimiter='\t')

    # create datetime variable
    datestr = '20'+df['YR'].map(str)+'-'+df['MO'].map(str)+'-'+df['DA'].map(str)
    timestr = df['HH'].map(str)+':'+df['MM'].map(str)+':'+df['SS'].map(str)+'.'+df['HH.1'].map(str)
    df['datetime'] = pd.to_datetime(datestr+' '+timestr,utc=True)

    hrow, = np.where(['Pings/Ens' in s for s in headerlines])
    PingsPerEns = int(headerlines[hrow.squeeze()].split('\t')[1])

    hrow, = np.where(['Time/Ping' in s for s in headerlines])
    TimePerPing = headerlines[hrow.squeeze()].split(' = ')[1]

    hrow, = np.where(['First Ensemble Date' in s for s in headerlines])
    FirstEnsDate = headerlines[hrow.squeeze()].split(' = ')[1]

    hrow, = np.where(['First Ensemble Time' in s for s in headerlines])
    FirstEnsTime = headerlines[hrow.squeeze()].split(' = ')[1]

    hrow, = np.where(['Ensemble Interval' in s for s in headerlines])
    EnsInterval = float(headerlines[hrow.squeeze()].split(' = ')[1])

    hrow, = np.where(['1st Bin Range' in s for s in headerlines])
    BlankDist = float(headerlines[hrow.squeeze()].split(' = ')[1])

    hrow, = np.where(['Bin Size' in s for s in headerlines])
    BinSize = float(headerlines[hrow.squeeze()].split('\t')[1])

    binnumbers = np.unique(headerlines[14].split('\t'))[1:]
    nbins = np.max([int(binnum) for binnum in binnumbers])

    BinHeight = np.arange(BlankDist,BlankDist+nbins*BinSize,BinSize)

    ntime = len(df['datetime'])
    nbeams = 4

    ds = xr.Dataset(coords={'time': df['datetime'].values})

    all_singlevars = ["Pit","Rol","Hea","Tem","Dep","Ori","BIT","Bat"]

    for var in all_singlevars:
        ds[var] = ('time', df[var].values)

    all_binvars = ['Eas','Nor','Ver','Err']

    for var in all_binvars:
        vardata = np.nan*np.ones([ntime,nbins])
        for n in range(nbins):
            if n == 0:
                col = var
            else:
                col = var+'.'+str(n)
            vardata[:,n] = df[col].values
        ds[var] = (['time','bin'],vardata)

    all_beamvars = ['EA','PG','C']

    for var in all_beamvars:
        vardata = np.nan*np.ones([ntime,nbins,nbeams])
        for beam in range(1,nbeams+1):
            for n in range(nbins):
                if n == 0:
                    col = var+str(beam)
                else:
                    col = var+str(beam)+'.'+str(n)
                vardata[:,n,beam-1] = df[col].values
            ds[var] = (['time','bin','beam'],vardata)
            
    # height of bins relative to instrument
    ds['binheight'] = ('bin',BinHeight)
    
    # Convert to units of m/s
    ds['Eas'] = ds['Eas']/1000
    ds['Nor'] = ds['Nor']/1000
    
    ds['Eas'].attrs['units'] = 'm/s'
    ds['Nor'].attrs['units'] = 'm/s'
    
    if time_span is not None:
        t1 = time_span[0]
        t2 = time_span[1]
        ii, = np.where((ds['time'] >= t1) & (ds['time'] <= t2))
        ds = ds.isel(time=ii)
    
    # add metadata
    ds.attrs['PingsPerEns'] = PingsPerEns
    ds.attrs['TimePerPing'] = TimePerPing    
    ds.attrs['First Ensemble Date'] = FirstEnsDate 
    ds.attrs['First Ensemble Time'] = FirstEnsTime
    ds.attrs['Ensemble Interval'] = EnsInterval
    ds.attrs['1st Bin Range'] = BlankDist
    ds.attrs['Bin Size'] = BinSize
    ds.attrs['RDI binary file'] = headerlines[1]
    ds.attrs['Instrument'] = headerlines[2]
            
    return ds