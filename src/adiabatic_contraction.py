import numpy as np
import sys
from scipy import interpolate
from parameters import *
from profiles import *
from astropy import units, constants


G = constants.G
G = G.to(units.kpc**3/units.Gyr**2.0/units.Msun)

# reading contra output file
contra_out = np.loadtxt(contra_output)
r_f = contra_out[:,1] * Rvir_host
rho_f = contra_out[:,5] * M_host / Rvir_host**3.0
pot_i = -167300.72
dpot_i = 7807.67

f1 = interpolate.interp1d(r_f, rho_f)
r_grid = np.linspace(1, Rvir_host, 1e6)
dr = r_grid[1] - r_grid[0]
rho_grid = f1(r_grid)
pot_ac = np.zeros(len(rho_grid))
dpot = np.zeros(len(rho_grid))
pot_ac[0] = pot_i
dpot[0] = dpot_i

for i in range(1,len(rho_grid)):
    pot_ac[i] = pot_ac[i-1] + dr*dpot[i-1]
    dpot[i] = dpot[i-1] + dr*(4.*np.pi*G.value*rho_grid[i] - 2./r_grid[i]*dpot[i-1])


    

def rho_ac(r):
    return f1(r)

def acc_ac(x, y, z):
    """
    Function that computes the acceleration of the satellite in a Halo 
    with Adiabatic Contraction (AC).

    Inputs:
    -------
    x, y, z
    Coordinates of the satellite position.

    Returns:
    --------
    the acceleration vector is cartesian coordinates.
    """
    r_eval = np.sqrt(x**2.+y**2.+z**2.)
    if r_eval > Rvir_host:
        print 'Warning: computing the acceleration with AC outside Rvir'
        sys.exit()

    f2 = interpolate.interp1d(r_grid, dpot)
    d_pot2 = f2(r_eval)

    theta = np.arccos(z/r_eval)
    phi = np.arctan2(y, x)
    ax = -d_pot2*np.sin(theta)*np.cos(phi)
    ay = -d_pot2*np.sin(theta)*np.sin(phi)
    az = -d_pot2*np.cos(theta)
    #print az, a_NFWnRvir(9.86, x, y, z, M_host, Rvir_host)[2], r_eval
    return ax, ay, az

"""
if __name__ == "__main__":

    g = constants.g
    g = g.to(units.kpc**3/units.gyr**2.0/units.msun)

    # reading contra output file
    contra_out = np.loadtxt(contra_output)
    r_f = contra_out[:,1] * rvir_host
    rho_f = contra_out[:,5] * m_host / rvir_host**3.0
    pot_i = -167300.72
    dpot_i = 7807.67

    f1 = interpolate.interp1d(r_f, rho_f)
    r_grid = np.linspace(1, rvir_host, 1e6)
    dr = r_grid[1] - r_grid[0]
    rho_grid = f1(r_grid)
    pot_ac = np.zeros(len(rho_grid))
    dpot = np.zeros(len(rho_grid))
    pot_ac[0] = pot_i
    dpot[0] = dpot_i

    for i in range(1,len(rho_grid)):
        pot_ac[i] = pot_ac[i-1] + dr*dpot[i-1]
        dpot[i] = dpot[i-1] + dr*(4.*np.pi*g.value*rho_grid[i] - 2./r_grid[i]*dpot[i-1])
"""

