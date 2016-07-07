import numpy as np
from scipy import interpolate

# Reading contra output file
contra_output = np.loadtxt(contra_output)
r_f = contra_output[:,1]
rho_f = contra_output[:,5]

def acceleration_components(r, pot):

def poisson_solver(pot_i, dpot_i, dr, r, rho):
    pot_ac = np.zeros(len(rho))
    dpot = np.zeros(len(rho))
    pot_ac[0] = pot_i
    dpot[0] = dpot_i
    for i in range(1,len(r)):
        pot_ac[i] = pot[i-1] + dr*dpot[i-1]
        dpot[i] = dpot[i-1] + dr*(4*np.pi*G*rho[i] - 2./r[i]*dpot[i-1])
    # truncate the halo at rvir
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
    f1 = interpolate.interp1d(r_f, rho_f)
    r_grid = np.linspace(0.001, 1, 10000)
    rho_grid = f1(r_grid)
    pot_grid = poisson_solve(rho_grid)
    f2 = interpolate.interp1d(r_grid, pot_grid)
    return ax, ay, az
    
