import numpy as np
from scipy import interpolate

# Reading contra output file
contra_output = np.loadtxt(contra_output)
r_f = contra_output[:,1]
rho_f = contra_output[:,5]

def integral():


def potential_ac(r):
    f = interpolate.interp1d(r_f, rho_f)
    
