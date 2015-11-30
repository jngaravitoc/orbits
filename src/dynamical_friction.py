import numpy as np
from scipy.special import erf
from astropy import units, constants
from profiles import *
from parameters import *

def coulomb_log(r):
    bmax = r # position of test particle at a given time
    # k = softening length if the satellite is modeled with a plummer
    # profile. See http://adsabs.harvard.edu/abs/2007ApJ...668..949B
    k = 3  #kpc
    bmin = 1.4 * k
    L = bmax / bmin
    # alpha is for make the dynamical friction more realistic.
    CL = alpha_df * np.log(L)
    return CL

def sigma(c, r, M, Rv):
    M = M * units.Msun
    Rv = Rv * units.kpc
    vvir = np.sqrt(G * M / Rv)
    g = np.log(1+c) - (c /(1+c))
    vmax = np.sqrt(0.216 * vvir**2 * c / g)
    x = r / Rv.value * c
    sigma = vmax * 1.4393 * x **(0.354) / (1 + 1.1756*x**0.725)
    sigma = sigma.to(units.kpc / units.Gyr)
    return sigma.value

def df(x, y, z, vx, vy, vz, M1, M2, Rv, c):
    # M2 would be the galaxy feeling the dynamical friction due to M1
    # Rv, c, (x, y, z) and (vx, vy, vz) are for the M2 galaxy
    # Coordinates
    r = np.sqrt(x**2 + y**2 + z**2)
    # Velocities
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    # Density of the NFW at a given r
    rho = dens_NFWnRvir(c, x, y, z, M1, Rv)
    # Mass of the satellite
    # Computing the dynamical friction
    factor = - 4 * np.pi * G**2
    factor = factor.to(units.kpc**6 / units.Msun**2 / units.Gyr**4)
    factor = factor.value
    Coulomb =  coulomb_log(r)
    s = sigma(c, r, M1, Rv)
    X = v / ( np.sqrt(2) * s)
    # Main equation
    a_dfx = (factor * M2 * rho * Coulomb * (erf(X) - 2.0 * X / (np.sqrt(np.pi)) * np.exp(-X**2.0)) * vx) / v**3.0
    a_dfy = (factor * M2 * rho * Coulomb * (erf(X) - 2.0 * X / (np.sqrt(np.pi)) * np.exp(-X**2.0)) * vy) / v**3.0
    a_dfz = (factor * M2 * rho * Coulomb * (erf(X) - 2.0 * X / (np.sqrt(np.pi)) * np.exp(-X**2.0)) * vz) / v**3.0
    # Units
    # kpc/Gyr2
    #print 'here inside the df function'
    #print a_dfx
    return a_dfx, a_dfy, a_dfz
