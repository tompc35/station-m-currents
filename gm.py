# Copyright (C) 2016 Joern Callies
#
# This file is part of GM81.
#
# GM81 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GM81 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GM81. If not, see <http://www.gnu.org/licenses/>.

# This module implemets the empirical spectrum of internal waves developed by
# Garrett and Munk, in the incarnation presented in Munk's chapter in Evolution
# of Physical Oceanography, which can be downloaded here:
# http://ocw.mit.edu/resources/res-12-000-evolution-of-physical-oceanography-spring-2007/part-2/wunsch_chapter9.pdf
# The variable names follow Munk's notation.

import numpy as np

# energy parameter
E = 6.3e-5

# j_*
js = 3

# sum_1^infty (j^2+j_*^2)^-1
jsum = (np.pi*js/np.tanh(np.pi*js)-1)/(2*js**2)

# gravitational acceleration
g = 9.81

def omg_k_j(k, j, f, N0, b):
    # frequency omega as a function of hor. wavenumber k and mode number j
    return np.sqrt((N0**2*k**2+f**2*(np.pi*j/b)**2)/(k**2+(np.pi*j/b)**2))

def k_omg_j(omg, j, f, N0, b):
    # hor. wavenumber as a function of frequency omg and mode number j
    return np.sqrt((omg**2-f**2)/(N0**2-omg**2))*np.pi*j/b

def B(omg, f):
    # Munk's B(omg) describing the frequency distribution
    return 2/np.pi*f/omg/np.sqrt(omg**2-f**2)

def H(j):
    # Munk's H(j) describing the mode distribution
    return 1./(j**2+js**2)/jsum

def E_omg_j(omg, j, f):
    # Munk's E(omg,j)
    return B(omg, f)*H(j)*E

def E_k_j(k, j, f, N, N0, b):
    # Munk's E(omg,j) transformed into hor. wavenumber space:
    # E(k,j) = E(omg,j) domg/dk. The transformation is achieved using the
    # dispersion relation (9.23a) in Munk (1981).
    omg = omg_k_j(k, j, f, N0, b)
    domgdk = (N0**2-omg**2)/omg*k/(k**2+(np.pi*j/b)**2)
    return E_omg_j(omg, j, f)*domgdk

def P_k_j(k, j, f, N, N0, b):
    # potential energy spectrum (N^2 times displacement spectrum) as a function
    # of hor. wavenumber k and mode number j
    omg = omg_k_j(k, j, f, N0, b)
    return b**2*N0*N*(omg**2-f**2)/omg**2*E_k_j(k, j, f, N, N0, b)

def K_k_j(k, j, f, N, N0, b):
    # kinetic energy spectrum as a function of hor. wavenumber k and mode
    # number j
    omg = omg_k_j(k, j, f, N0, b)
    return b**2*N0*N*(omg**2+f**2)/omg**2*E_k_j(k, j, f, N, N0, b)

def eta_k_j(k, j, f, N, N0, b):
    # SSH spectrum as a function of hor. wavenumber k and mode number j
    omg = omg_k_j(k, j, f, N0, b)
    return (omg**2-f**2)**2/(f**2*(omg**2+f**2))*K_k_j(k, j, f, N, N0, b)/k**2*f**2/g**2

def P_omg_j(omg, j, f, N, N0, b):
    # potential energy spectrum (N^2 times displacement spectrum) as a function
    # of frequency omg and mode number j
    return b**2*N0*N*(omg**2-f**2)/omg**2*E_omg_j(omg, j, f)

def K_omg_j(omg, j, f, N, N0, b):
    # kinetic energy spectrum as a function of frequency omg and mode number j
    return b**2*N0*N*(omg**2+f**2)/omg**2*E_omg_j(omg, j, f)

def eta_omg_j(omg, j, f, N, N0, b):
    # SSH spectrum as a function of frequency omg and mode number j
    k = k_omg_j(omg, j, f, N0, b)
    return (omg**2-f**2)**2/(f**2*(omg**2+f**2))*K_omg_j(omg, j, f, N, N0, b)/k**2*f**2/g**2

def sqrt_trapz(kh, S):
    # integrate S/sqrt(kh^2-k^2) over all kh, approximating S as piecewise
    # linear but then performing the integration exactly
    a = kh[:-1]
    b = kh[1:]
    A = S[:-1]
    B = S[1:]
    k = kh[0]
    return np.sum(((A-B)*(np.sqrt(a**2-k**2)-np.sqrt(b**2-k**2))+(a*B-b*A)*np.log((a+np.sqrt(a**2-k**2))/(b+np.sqrt(b**2-k**2))))/(b-a))

def calc_1d(k, S):
    # calculate 1D wavenumber spectrum from 2D isotropic wavenumber spectrum:
    # S1d = 2/pi int_k^infty S2d/sqrt(kh^2-k^2) dkh
    # (The normalization is such that int_0^infty S1d dk = int_0^infty S2d dkh.)
    S1d = np.empty(k.size)
    for i in range(k.size):
        S1d[i] = 2/np.pi*sqrt_trapz(k[i:], S[i:])
    return S1d
