import numpy as np

def cor_beta(lat):
    '''
    Calculate the meridional gradient of the Coriolis frequency, beta = df/dy
    
    INPUT
    lat: latitude (degrees)
    
    RETURNS
    beta (m^-1 s^-1)
    '''
    omega = 7.2921e-5
    earth_radius = 6.371e6
    return 2 * omega * np.cos(lat*np.pi/180) / earth_radius

def bottomstress(u,v,z,zo=1e-3,rho=1025):
    '''
    Calculate bottom stress using logarithmic boundary layer formulation.

    INPUT
    u,v: components of near-bottom velocity (m/s)
    z: height off bottom (m, positive)
    zo: roughness (m, default 1e-3)
    rho: density (kg/m^3, default 1025)
    
    RETURNS
    taux: x-component of stress (N m^-2)
    tauy: y-component of stress (N m^-2)
    '''

    kappa = 0.41 # Von Karman's constant

    # speed and orientation
    umag = np.abs(u+1j*v)
    theta = np.angle(u+1j*v)

    # friction velocity
    ustar = umag/((1/kappa)*np.log(z/zo))
    taumag = rho*ustar**2

    # vector components
    taux = -taumag*np.cos(theta)
    tauy = -taumag*np.sin(theta)

    return taux,tauy
