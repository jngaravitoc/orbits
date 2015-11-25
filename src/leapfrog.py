import numpy as np
from astropy import units, constants
from acceleration import *

def leapfrog(x_ic, y_ic, z_ic, vx_ic, vy_ic, vz_ic, M_halo, M_disk, M_bulge, M_sat, Rvir):

    n_points = 2500
    h = 0.001
    # Creating the arrays to collect the data in each step of the integration
    # the imput units should be in Kpc for positions and km/s for velocities

    # Coordinates tranformation

    vx_ic = vx_ic * units.km / units.s
    vy_ic = vy_ic * units.km / units.s
    vz_ic = vz_ic * units.km / units.s
    vx_ic = vx_ic.to(units.kpc / units.Gyr)
    vy_ic = vy_ic.to(units.kpc / units.Gyr)
    vz_ic = vz_ic.to(units.kpc / units.Gyr)
    vx_ic = vx_ic.value
    vy_ic = vy_ic.value
    vz_ic = vz_ic.value

    t = np.zeros(n_points)
    x = np.zeros(n_points) #+ x_ic
    y = np.zeros(n_points) #+ y_ic
    z = np.zeros(n_points) #+ z_ic

    x_mw = np.zeros(n_points)
    y_mw = np.zeros(n_points)
    z_mw = np.zeros(n_points)


    vx = np.zeros(n_points)
    vy = np.zeros(n_points)
    vz = np.zeros(n_points)

    vx_mw  = np.zeros(n_points)
    vy_mw  = np.zeros(n_points)
    vz_mw  = np.zeros(n_points)

    vx_ic_MW = 0 * units.km / units.s
    vy_ic_MW = 0 * units.km / units.s
    vz_ic_MW = 0 * units.km / units.s

    vx_ic_MW = vx_ic_MW.to(units.kpc / units.Gyr)
    vy_ic_MW = vy_ic_MW.to(units.kpc / units.Gyr)
    vz_ic_MW = vz_ic_MW.to(units.kpc / units.Gyr)

    ax = np.zeros(n_points)
    ay = np.zeros(n_points)
    az = np.zeros(n_points)

    ax_mw = np.zeros(n_points)
    ay_mw = np.zeros(n_points)
    az_mw = np.zeros(n_points)

    t[0] = 0

    # This initial conditions come form MW.py, the units are Kpc and Gyr
    x[0] = x_ic
    y[0] = y_ic
    z[0] = z_ic

    x_mw[0] = 0
    y_mw[0] = 0
    z_mw[0] = 0

    vx[0] = vx_ic
    vy[0] = vy_ic
    vz[0] = vz_ic

    vx_mw[0] = vx_ic_MW.value
    vy_mw[0] = vy_ic_MW.value
    vz_mw[0] = vz_ic_MW.value

    ax[0] = acceleration((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]), (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]), M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
    ay[0] = acceleration((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]), (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]), M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
    az[0] = acceleration((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]), (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]), M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

    ax_mw[0] = acceleration_mw((x_mw[0]-x[0]),-(y[0] - y_mw[0]), -(z[0] - z_mw[0]), (vx_mw[0] - vx[0]), (vy_mw[0] - vy[0]), (vz_mw[0] - vz[0]), M_sat, M_halo)[0]
    ay_mw[0] = acceleration_mw((x_mw[0]-x[0]),-(y[0] - y_mw[0]), -(z[0] - z_mw[0]), (vx_mw[0] - vx[0]), (vy_mw[0] - vy[0]), (vz_mw[0] - vz[0]), M_sat, M_halo)[1]
    az_mw[0] = acceleration_mw((x_mw[0]-x[0]),-(y[0] - y_mw[0]), -(z[0] - z_mw[0]), (vx_mw[0] - vx[0]), (vy_mw[0] - vy[0]), (vz_mw[0] - vz[0]), M_sat, M_halo)[2]


    # one half step
    # at the first step the MW is at 0, 0, 0 and its not moving

    t[1] = t[0] - h
    x[1] = x[0] - h * vx[0]
    y[1] = y[0] - h * vy[0]
    z[1] = z[0] - h * vz[0]

    vx[1] = vx[0] - h * ax[0]
    vy[1] = vy[0] - h * ay[0]
    vz[1] = vz[0] - h * az[0]
    if (mw=='free'):
        x_mw[1] = x_mw[0] - h * vx_mw[0]
        y_mw[1] = y_mw[0] - h * vy_mw[0]
        z_mw[1] = z_mw[0] - h * vz_mw[0]

        vx_mw[1] = vx_mw[0] - h * ax_mw[0]
        vy_mw[1] = vy_mw[0] - h * ay_mw[0]
        vz_mw[1] = vz_mw[0] - h * az_mw[0]

        ax_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), (vx_mw[1] - vx[1]), (vy_mw[1] - vy[1]), (vz_mw[1] - vz[1]), M_sat, M_halo)[0]
        ay_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), (vx_mw[1] - vx[1]), (vy_mw[1] - vy[1]), (vz_mw[1] - vz[1]), M_sat, M_halo)[1]
        az_mw[1] = acceleration_mw(-(x[1] - x_mw[1]), -(y[1] - y_mw[1]), -(z[1] - z_mw[1]), (vx_mw[1] - vx[1]), (vy_mw[1] - vy[1]), (vz_mw[1] - vz[1]), M_sat, M_halo)[2]



    ax[1] = acceleration(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
    ay[1] = acceleration(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
    az[1] = acceleration(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

    for i in range(2,n_points):
        t[i] = t[i-1] - h

        x[i] = x[i-2] - 2 * h * vx[i-1]
        y[i] = y[i-2] - 2 * h * vy[i-1]
        z[i] = z[i-2] - 2 * h * vz[i-1]

        vx[i] = vx[i-2] - 2 * h * ax[i-1]
        vy[i] = vy[i-2] - 2 * h * ay[i-1]
        vz[i] = vz[i-2] - 2 * h * az[i-1]

        if (mw=='free'):
            x_mw[i] = x_mw[i-2] - 2 * h * vx_mw[i-1]
            y_mw[i] = y_mw[i-2] - 2 * h * vy_mw[i-1]
            z_mw[i] = z_mw[i-2] - 2 * h * vz_mw[i-1]

            vx_mw[i] = vx_mw[i-2] - 2 * h * ax_mw[i-1]
            vy_mw[i] = vy_mw[i-2] - 2 * h * ay_mw[i-1]
            vz_mw[i] = vz_mw[i-2] - 2 * h * az_mw[i-1]


            ax_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), (vx_mw[i] - vx[i]), (vy_mw[i] - vy[i]), (vz_mw[i] - vz[i]), M_sat, M_halo)[0]
            ay_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), (vx_mw[i] - vx[i]), (vy_mw[i] - vy[i]), (vz_mw[i] - vz[i]), M_sat, M_halo)[1]
            az_mw[i] = acceleration_mw((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), (vx_mw[i] - vx[i]), (vy_mw[i] - vy[i]), (vz_mw[i] - vz[i]), M_sat, M_halo)[2]


            ax[i] = acceleration(x[i-1]-x_mw[i-1], y[i-1]-y_mw[i-1], z[i-1]-z_mw[i-1], vx[i-1]-vx_mw[i-1], vy[i-1]-vy_mw[i-1], vz[i-1]-vz_mw[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[0]
            ay[i] = acceleration(x[i-1]-x_mw[i-1], y[i-1]-y_mw[i-1], z[i-1]-z_mw[i-1], vx[i-1]-vx_mw[i-1], vy[i-1]-vy_mw[i-1], vz[i-1]-vz_mw[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[1]
            az[i] = acceleration(x[i-1]-x_mw[i-1], y[i-1]-y_mw[i-1], z[i-1]-z_mw[i-1], vx[i-1]-vx_mw[i-1], vy[i-1]-vy_mw[i-1], vz[i-1]-vz_mw[i-1], M_halo, M_disk, M_bulge, M_sat, Rvir)[2]

        return t, x, y, z, x_mw, y_mw, z_mw, vx, vy, vz, vx_mw, vy_mw, vz_mw, ax, ay, az, ax_mw, ay_mw, az_mw
