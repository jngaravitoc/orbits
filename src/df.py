import numpy as np
from scipy import erf
from astropy import units, constants


def coulomb_log(r):
    bmax = r # position of test particle at a time t
    k = 3 * units.kpc # kpc
    bmin = 1.4 * k # k is the softening length if the LMC were modeled using a plummer profile . See Besla07
    L = bmax / bmin
    CL = alpha * np.log(L)
    return CL

def sigma(r_s, r, M_halo, Rv):
    M_halo = M_halo * units.Msun
    Rvir = Rv * units.kpc 
    vvir = np.sqrt(G * M_halo / Rvir) 
    c = Rvir.value / r_s
    g = np.log(1+c) - (c /(1+c))
    vmax = np.sqrt(0.216 * vvir**2 * c / g)
    x = r.value / rs
    sigma = vmax * 1.4393 * x **(0.354) / (1 + 1.1756*x**0.725)
    sigma = sigma.to(units.kpc / units.Gyr)
    return sigma

def dynamical_friction_sis(x, y, z, vx, vy, vz, M_halo, M_sat, Rvir, r_s):
    # Coordinates
    x = x * units.kpc
    y = y * units.kpc
    z = z * units.kpc
    r = np.sqrt(x**2 + y**2 + z**2)
    # Velocities
    vx = vx * units.kpc / units.Gyr
    vy = vy * units.kpc / units.Gyr
    vz = vz * units.kpc / units.Gyr
    v = np.sqrt(vx**2 + vy**2 + vz**2)
    # Density of the NFW at a given r
    Mhalo = M_halo
    c = Rvir / r_s
    rho = dens_NFWnRvir(c, x.value, y.value, z.value, Mhalo, Rvir) 
    # Mass of the satellite
    M_sat = M_sat * units.Msun
    # Computing the dynamical friction
    #G = constants.G
    #G = G.to(units.kpc**3 / units.Msun / units.Gyr**2)
    factor = - 4 * np.pi * G**2
    Coulomb =  coulomb_log(r) #**************
    s = sigma(r_s, r, M_halo, Rvir) 
    X = v / ( np.sqrt(2) * s ) 
    #print v, s, X
    # Main equation
    a_dfx = (factor * M_sat * rho * Coulomb  * (  erf(X) - 2.0*X/(np.sqrt(np.pi)) * np.exp(-X**2.0)  ) * vx) / v**3.0
    a_dfy = (factor * M_sat * rho * Coulomb  * (  erf(X) - 2.0*X/(np.sqrt(np.pi)) * np.exp(-X**2.0)  ) * vy) / v**3.0
    a_dfz = (factor * M_sat * rho * Coulomb  * (  erf(X) - 2.0*X/(np.sqrt(np.pi)) * np.exp(-X**2.0)  ) * vz) / v**3.0 
    # Transforming to the right units
    a_dfx = a_dfx.to(units.kpc / units.Gyr**2) 
    a_dfy = a_dfy.to(units.kpc / units.Gyr**2)
    a_dfz = a_dfz.to(units.kpc / units.Gyr**2)
    return a_dfx.value, a_dfy.value, a_dfz.value
