import numpy as np

def cor_beta(lat):
    omega = 7.2921e-5
    earth_radius = 6.371e6
    return 2 * omega * np.cos(lat*np.pi/180) / earth_radius