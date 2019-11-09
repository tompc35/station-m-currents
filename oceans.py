import numpy as np

def cor_beta(lat):
    omega = 7.2921e-5
    earth_radius = 6.371e6
    return 2 * omega * np.cos(lat*np.pi/180) / earth_radius

def bottomstress(u,v,z,zo=1e-3,rho=1025):
    '''
    Calculate bottom stress using logarithmic boundary layer formulation.

    u,v: components of near-bottom velocity (m/s)
    z: height off bottom (m, positive)
    zo: roughness (m)
    rho: density (kg/m^3, default)
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
