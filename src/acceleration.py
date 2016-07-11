"""
This code computes the acceleration for the satellite and the host.

"""

import numpy as np
from astropy import units, constants
from profiles import *
from parameters import *
from dynamical_friction import *
from adiabatic_contraction import *


#Function that computes the satellite acceleration
def acc_sat(x, y, z, vx, vy, vz):
    r = np.sqrt(x**2 + y**2 + z**2)
    # Acceleration by the DM halo profile
    if (r<=Rvir_host):
        if (ac == 1):
            ahalo = acc_ac(x, y, z)
        else:
            if (Host_model == 0):
                ahalo = a_NFWnRvir(c_host, x, y, z, M_host, Rvir_host)
            elif (Host_model == 1):
                ahalo = a_hernquist(rs_host, x, y, z, M_host)

        adisk = a_mn(a_disk, b_disk, x, y, z, M_disk)
        abulge = a_hernquist(rh, x, y, z, M_bulge)

        ax = ahalo[0] + adisk[0] + abulge[0]
        ay = ahalo[1] + adisk[1] + abulge[1]
        az = ahalo[2] + adisk[2] + abulge[2]

        # Truncating the halo at r_vir:
        # Dynamical Friction inside the r_vir
        a_dfx, a_dfy, a_dfz = df(x, y, z, vx, vy, vz, M_host, M_sat, \
                              Rvir_host, c_host, alpha_df_sat)
        Ax = ax + a_dfx
        Ay = ay + a_dfy
        Az = az + a_dfz
    # Point like acceleration beyond r_vir
    else:
        Mtot = (M_host + M_disk + M_bulge) * units.Msun
        Ax = - G * Mtot * x * units.kpc / (r*units.kpc)**3
        Ay = - G * Mtot * y * units.kpc / (r*units.kpc)**3
        Az = - G * Mtot * z * units.kpc / (r*units.kpc)**3
        Ax = Ax.to(units.kpc / units.Gyr**2)
        Ay = Ay.to(units.kpc / units.Gyr**2)
        Az = Az.to(units.kpc / units.Gyr**2)
        Ax = Ax.value
        Ay = Ay.value
        Az = Az.value
    return Ax, Ay, Az

# Function that computes the host acceleration
def acc_host(x, y, z, vx, vy, vz):
    #Host halo acceleration due to the Satellite DM halo profile.
    if (Sat_model == 2):
        A_host = a_NFWnRvir(c_sat, x, y, z, M_sat, Rvir_sat)

    elif (Sat_model == 1):
        A_host = a_hernquist(rs_sat, x, y, z, M_sat)

    elif (Sat_model == 0):
        A_host = a_plummer(rs_sat, x, y, z, M_sat)

    Ax = A_host[0]
    Ay = A_host[1]
    Az = A_host[2]

    # Host dynamical friction:
    if (Host_df==1):
        D = np.sqrt(x**2 + y**2 + z**2)
        R_mass = Rvir_host - (Rvir_sat - D)
        # Mass fraction of the host galaxy inside the satellite.
        if (Host_model == 0):
            M_frac = mass_NFWnRvir(c_host, R_mass, 0, 0, M_host, Rvir_host)
        elif (Host_model ==1):
            M_frac = mass_hernquist(rs_host, R_mass,M_host)
        a_dfx, a_dfy, a_dfz = df(x, y, z, vx, vy, vz, M_sat, M_frac,\
                              Rvir_sat, c_sat, alpha_df_host)
        Ax = Ax + a_dfx
        Ay = Ay + a_dfy
        Az = Az + a_dfz

    return Ax, Ay, Az
