import numpy as np
from scipy.signal import get_window

def lombscargle(t,x,ofac=4,hifac=1,t0=None,return_onesided=True,return_zero=False, window='boxcar',scaling='classical'):
    '''
    Compute the discrete Fourier transform and periodogram for unevenly-spaced 
    data using the Lomb-Scargle periodogram. Follows methods outlined in Scargle
    (1989). 
    
    INPUTS
    
    t - array of numerical time values (length N)
    x - array of data values, may be complex (length N)
    
    RETURNS
    
    f - array of frequencies
    ftx - array of complex coefficients
          discrete Fourier transform of x
    px - periodogram of ftx, proportional to |ftx|**2
         with the default "classical" scaling used in Scargle (1989), 
         px = (1/N)*|ftx|**2
    
    OPTIONAL PARAMETERS
    
    ofac - oversampling parameter 
           ratio of number of frequencies used to number of samples in x
           (default 4)        
    hifac - high frequency parameter
            ratio of highest frequency to pseudo-Nyquist frequency
            (default 1)  
    t0 - time origin
         reference point for phase calculation
         (default - None, first value in t array is used)
    return_onesided - boolean for returning a one-sided spectrum
        If True, return a one-sided spectrum for real data. 
        If False return a two-sided spectrum for real data.  
        Note that for complex data, a two-sided spectrum is always returned.
        (default - True)
    return_zero - boolean for evaluating zero frequency
        If True, include zero frequency. If False, do not include zero frequency.
        Uses expressions for the limit as frequency approaches zero, following
        Scargle (1989). 
        (default - False)
    window - String specifying desired window to use. See `scipy.signal.get_window` for 
        a list of windows and required parameters.
    scaling - Selects between computing the classical periodogram used by Scargle 
        ('classical') or the power spectral density ('density'). The classical 
        periodogram  has units of x**2. The power  spectral density has units of x**2/f. 
        The scaling determines how the periodogram px is calculated from the discrete 
        Fourier transform ftx:
        'classical': px = (1/N)*|ftx|**2, where N is the number of samples
        'density': px = (deltat/N)*|ftx|**2, where deltat is the average time step
      
    REFERENCE
    
    Scargle, J.D. (1989) Studies in astronomical time series analysis III: Fourier transforms,
    autocorrelation functions, and cross-correlation functions of unevenly spaced data. The
    Astrophysical Journal, 343, 874-887
    '''

    i = 1j # square root of -1
    N = len(x) # number of samples    
    
    wts = get_window(window,N)
    wts = N*wts/np.sum(wts) # make sum of weights equal to N
    
    x = x*wts  # apply window

    intm = np.mean(np.diff(t))

    flo = ((intm)**-1)/(len(x)*ofac)  # lowest freq
    fhi = hifac*(2*intm)**-1          # highest freq

    f = np.arange(flo,fhi+flo,flo)
    
    if return_zero == True:
        f = np.append(0,f)
        
    # if complex, evaluate two-sided spectrum regardless of user choice
    if np.any(np.iscomplex(x)):
        return_onesided = False
    
    # two-sided spectrum
    if return_onesided == False:
        if return_zero == True:
            f = np.append(-f[1:][::-1],f)
        else:
            f = np.append(-f[::-1],f)

    # time origin (reference point for phase calculation)
    if t0 is None:
        t0 = t[0] 

    # initialize DFT as array of complex numbers
    ftx = np.nan*np.ones(len(f)) + i*np.nan*np.ones(len(f))

    for k,fk in enumerate(f):
        wrun = 2*np.pi*fk # angular frequency    
        
        if fk == 0:
            # use well-defined limit as frequency approaches zero
            tau = np.sum(t)/N
            ftx[k] = np.sum(x)/np.sqrt(N)
            
        else:
            Fo = ((N/2)**0.5)*np.exp(-i*wrun*t0) 

            tau = np.arctan2(np.sum(np.sin(2*wrun*t)),np.sum(np.cos(2*wrun*t)))/(2*wrun)
            tprime = t - tau

            A = np.sum(np.cos(wrun*tprime)**2)**-0.5
            B = np.sum(np.sin(wrun*tprime)**2)**-0.5

            # Note apparent typo in Scargle (1989), which has a plus 
            # sign (+) instead of a minus sign below. This only makes a 
            # difference in the periodogram if the input values in x are complex.
            
            ftx[k] = Fo*np.sum(A*x*np.cos(wrun*tprime) - i*B*x*np.sin(wrun*tprime))

    if scaling == 'classical':
        px = np.abs(ftx)**2/N
    elif scaling == 'density':
        px = np.abs(ftx)**2/N*intm
    else:
        raise ValueError('Scaling argument not understood. Acceptable options are classical or density')
    
    return f, ftx, px

# def veccor(u1,v1,u2,v2):
#     '''
    
#     RETURNS:
    
#     '''
    
    
        