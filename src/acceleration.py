import numpy as np
from astropy import units, constants
from profiles import *

def acc_sat(x, y, z, vx, vy, vz, M_halo, M_disk, M_bulge, M_sat, Rvir):
    c = Rvir / rs
    ahalo = a_NFWnRvir(c, x, y, z, M_halo, Rvir)
    adisk = a_mn(ra, rb, x, y, z, M_disk)
    abulge = a_hernquist(0.7, x, y, z, M_bulge)
    ax = ahalo[0] + adisk[0] + abulge[0]
    ay = ahalo[1] + adisk[1] + abulge[1]
    az = ahalo[2] + adisk[2] + abulge[2]
    ax = ax.to(units.kpc/units.Gyr**2)  
    ay = ay.to(units.kpc/units.Gyr**2) 
    az = az.to(units.kpc/units.Gyr**2)
    r = np.sqrt(x**2 + y**2 + z**2)
    # Truncating the halo at the virial radius

    if (r <= Rvir):
        a_dfx, a_dfy, a_dfz = dynamical_friction_sis(x, y, z, vx, vy, vz, M_halo, M_sat, Rvir, rs)
        Ax = ax.value + a_dfx
        Ay = ay.value + a_dfy
        Az = az.value + a_dfz
    else:
        Mtot = (M_halo) * units.Msun
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

def acc_host(x, y, z, vx, vy, vz, M_sat, M_halo):
    M_halo = M_halo * units.Msun
    M_sat = M_sat * units.Msun
    r = np.sqrt(x**2 + y**2 + z**2)
    Ax =  a_NFWnRvir(c_sat, x, y, z, M_sat.value, Rvir_sat)[0]
    Ay =  a_NFWnRvir(c_sat, x, y, z, M_sat.value, Rvir_sat)[1]
    Az =  a_NFWnRvir(c_sat, x, y, z, M_sat.value, Rvir_sat)[2]

    Ax = Ax.to(units.kpc / units.Gyr**2)
    Ay = Ay.to(units.kpc / units.Gyr**2)
    Az = Az.to(units.kpc / units.Gyr**2)

    if (r <= Rvir_sat):
        a_dfx, a_dfy, a_dfz = dynamical_friction_sis(x, y, z, vx, vy, vz, M_sat.value, M_halo.value, Rvir_sat, rs_sat)
        Ax = Ax.value + a_dfx
        Ay = Ay.value + a_dfy
        Az = Az.value + a_dfz
    return Ax, Ay, Az
