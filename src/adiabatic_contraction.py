import numpy as np
import sys
from scipy import interpolate
from parameters import *
from profiles import *

# Reading contra output file
contra_out = np.loadtxt(contra_output)

r_f = contra_out[:,1] * 3.*Rvir_host
rho_f = contra_out[:,5] * M_host / Rvir_host**3.0

pot_i = 0.
dpot_i = 0.

def poisson_solver(pot_i, dpot_i, dr, r, rho):
    pot_ac = np.zeros(len(rho))
    dpot = np.zeros(len(rho))
    pot_ac[0] = pot_i
    dpot[0] = dpot_i
    for i in range(1,len(rho)):
        pot_ac[i] = pot_ac[i-1] + dr*dpot[i-1]
        dpot[i] = dpot[i-1] + dr*(4*np.pi*G.value*rho[i] - 2./r[i]*dpot[i-1])
    return dpot



# Returns the acceleration taking into account adiabatic contraction.
def acc_ac(x, y, z):
    """
    Function that computes the acceleration of the satellite in a Halo 
    with Adiabatic Contraction (AC).
    
    Inputs:
    ------   
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

    f1 = interpolate.interp1d(r_f, rho_f)
    r_grid = np.linspace(1, Rvir_host, 1E6)
    dr = r_grid[1] - r_grid[0]
    rho_grid = f1(r_grid)

    d_pot = poisson_solver(pot_i, dpot_i, dr, r_grid, rho_grid)

    f2 = interpolate.interp1d(r_grid, d_pot)
    d_pot2 = f2(r_eval)

    theta = np.arccos(z/r_eval)
    phi = np.arctan2(y, x)
    ax = -d_pot2*np.sin(theta)*np.cos(phi)
    ay = -d_pot2*np.sin(theta)*np.cos(phi)
    az = -d_pot2*np.cos(theta)
    return ax, ay, az
