import numpy as np
from astropy import units, constants
from acceleration import *
from parameters import *

def leapfrog():
    # h is the time step
    h = 0.001 * direction
    n_points = int(time * 1000.0)

    t = np.zeros(n_points)
    x = np.zeros(n_points)
    y = np.zeros(n_points)
    z = np.zeros(n_points)

    x_mw = np.zeros(n_points)
    y_mw = np.zeros(n_points)
    z_mw = np.zeros(n_points)

    vx = np.zeros(n_points)
    vy = np.zeros(n_points)
    vz = np.zeros(n_points)

    vx_mw  = np.zeros(n_points)
    vy_mw  = np.zeros(n_points)
    vz_mw  = np.zeros(n_points)

    ax = np.zeros(n_points)
    ay = np.zeros(n_points)
    az = np.zeros(n_points)

    ax_mw  = np.zeros(n_points)
    ay_mw  = np.zeros(n_points)
    az_mw  = np.zeros(n_points)

    t[0] = 0
    x[0] = x_sat
    y[0] = y_sat
    z[0] = z_sat

    x_mw[0] = x_host
    y_mw[0] = y_host
    z_mw[0] = z_host

    vx[0] = vx_sat
    vy[0] = vy_sat
    vz[0] = vz_sat

    vx_mw[0] = vx_host
    vy_mw[0] = vy_host
    vz_mw[0] = vz_host

    ax[0] = acc_sat((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]),\
            (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]))[0]
    ay[0] = acc_sat((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]),\
            (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]))[1]
    az[0] = acc_sat((x[0]-x_mw[0]), (y[0]-y_mw[0]), (z[0]-z_mw[0]),\
            (vx[0]-vx_mw[0]), (vy[0]-vy_mw[0]), (vz[0]-vz_mw[0]))[2]

    if (sat2==True):
        ax2[0] = acc_sat((x2[0]-x_mw[0]), (y2[0]-y_mw[0]), (z2[0]-z_mw[0]),\
                (vx2[0]-vx_mw[0]), (vy2[0]-vy_mw[0]), (vz2[0]-vz_mw[0]))[0]
        ay2[0] = acc_sat2((x2[0]-x_mw[0]), (y2[0]-y_mw[0]), (z2[0]-z_mw[0]),\
                (vx2[0]-vx_mw[0]), (vy2[0]-vy_mw[0]), (vz2[0]-vz_mw[0]))[1]
        az2[0] = acc_sat2((x2[0]-x_mw[0]), (y2[0]-y_mw[0]), (z2[0]-z_mw[0]),\
                (vx2[0]-vx_mw[0]), (vy2[0]-vy_mw[0]), (vz2[0]-vz_mw[0]))[2]

    ax_mw[0] = acc_host((x_mw[0]-x[0]), (y_mw[0] - y[0]), (z_mw[0] -z[0]),\
               (vx_mw[0] - vx[0]), (vy_mw[0] - vy[0]), (vz_mw[0] - vz[0]))[0]
    ay_mw[0] = acc_host((x_mw[0]-x[0]), (y_mw[0] - y[0]), (z_mw[0] -z[0]),\
               (vx_mw[0] - vx[0]), (vy_mw[0] - vy[0]), (vz_mw[0] - vz[0]))[1]
    az_mw[0] = acc_host((x_mw[0]-x[0]), (y_mw[0] - y[0]), (z_mw[0] -z[0]),\
               (vx_mw[0] - vx[0]), (vy_mw[0] - vy[0]), (vz_mw[0] - vz[0]))[2]

    # half step
    # Here I assume the host galaxy starts at position (0, 0, 0) and then its
    # initial v[1] is (0, 0, 0)
    t[1] = t[0] - h
    x[1] = x[0] - h * vx[0]
    y[1] = y[0] - h * vy[0]
    z[1] = z[0] - h * vz[0]

    vx[1] = vx[0] - h * ax[0]
    vy[1] = vy[0] - h * ay[0]
    vz[1] = vz[0] - h * az[0]

    if (Host_move==1):
        x_mw[1] = x_mw[0] - h * vx_mw[0]
        y_mw[1] = y_mw[0] - h * vy_mw[0]
        z_mw[1] = z_mw[0] - h * vz_mw[0]

        vx_mw[1] = vx_mw[0] - h * ax_mw[0]
        vy_mw[1] = vy_mw[0] - h * ay_mw[0]
        vz_mw[1] = vz_mw[0] - h * az_mw[0]

        ax_mw[1] = acc_host((x_mw[1] - x[1]), (y_mw[1] - y[1]), (z_mw[1] - z[1]), (vx_mw[1] - vx[1]), (vy_mw[1] - vy[1]), (vz_mw[1] - vz[1]))[0]
        ay_mw[1] = acc_host((x_mw[1] - x[1]), (y_mw[1] - y[1]), (z_mw[1] - z[1]), (vx_mw[1] - vx[1]), (vy_mw[1] - vy[1]), (vz_mw[1] - vz[1]))[1]
        az_mw[1] = acc_host((x_mw[1] - x[1]), (y_mw[1] - y[1]), (z_mw[1] - z[1]), (vx_mw[1] - vx[1]), (vy_mw[1] - vy[1]), (vz_mw[1] - vz[1]))[2]

    ax[1] = acc_sat(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1])[0]
    ay[1] = acc_sat(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1])[1]
    az[1] = acc_sat(x[1]-x_mw[1], y[1]-y_mw[1], z[1]-z_mw[1], vx[1]-vx_mw[1], vy[1]-vy_mw[1], vz[1]-vz_mw[1])[2]



    if (sat2==True):
        ax2[1] = acc_sat((x2[1]-x_mw[1]), (y2[1]-y_mw[1]),(z2[1]-z_mw[1]),\
                (vx2[1]-vx_mw[1]), (vy2[1]-vy_mw[1]),(vz2[1]-vz_mw[1]))[0]
        ay2[1] = acc_sat2((x2[1]-x_mw[1]), (y2[1]-y_mw[1]), (z2[1]-z_mw[1]),\
                (vx2[1]-vx_mw[1]), (vy2[1]-vy_mw[1]), (vz2[1]-vz_mw[1]))[1]
        az2[1] = acc_sat2((x2[1]-x_mw[1]), (y2[1]-y_mw[1]), (z2[1]-z_mw[1]),\
                (vx2[1]-vx_mw[1]), (vy2[1]-vy_mw[1]), (vz2[1]-vz_mw[1]))[2]

    for i in range(2, len(x)):
        t[i] = t[i-1] - h
        x[i] = x[i-2] - 2 * h * vx[i-1]
        y[i] = y[i-2] - 2 * h * vy[i-1]
        z[i] = z[i-2] - 2 * h * vz[i-1]

        vx[i] = vx[i-2] - 2 * h * ax[i-1]
        vy[i] = vy[i-2] - 2 * h * ay[i-1]
        vz[i] = vz[i-2] - 2 * h * az[i-1]

        if (Host_move==1):
            x_mw[i] = x_mw[i-2] - 2 * h * vx_mw[i-1]
            y_mw[i] = y_mw[i-2] - 2 * h * vy_mw[i-1]
            z_mw[i] = z_mw[i-2] - 2 * h * vz_mw[i-1]

            vx_mw[i] = vx_mw[i-2] - 2 * h * ax_mw[i-1]
            vy_mw[i] = vy_mw[i-2] - 2 * h * ay_mw[i-1]
            vz_mw[i] = vz_mw[i-2] - 2 * h * az_mw[i-1]

            ax_mw[i] = acc_host((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), (vx_mw[i] - vx[i]), (vy_mw[i] - vy[i]), (vz_mw[i] - vz[i]))[0]
            ay_mw[i] = acc_host((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), (vx_mw[i] - vx[i]), (vy_mw[i] - vy[i]), (vz_mw[i] - vz[i]))[1]
            az_mw[i] = acc_host((x_mw[i]-x[i]), (y_mw[i] - y[i]), (z_mw[i]-z[i]), (vx_mw[i] - vx[i]), (vy_mw[i] - vy[i]), (vz_mw[i] - vz[i]))[2]

        ax[i] = acc_sat(x[i]-x_mw[i], y[i]-y_mw[i], z[i]-z_mw[i], vx[i]-vx_mw[i], vy[i]-vy_mw[i], vz[i]-vz_mw[i])[0]
        ay[i] = acc_sat(x[i]-x_mw[i], y[i]-y_mw[i], z[i]-z_mw[i], vx[i]-vx_mw[i], vy[i]-vy_mw[i], vz[i]-vz_mw[i])[1]
        az[i] = acc_sat(x[i]-x_mw[i], y[i]-y_mw[i], z[i]-z_mw[i], vx[i]-vx_mw[i], vy[i]-vy_mw[i], vz[i]-vz_mw[i])[2]

    return t, x, y, z, vx, vy, vz, x_mw, y_mw, z_mw, vx_mw, vy_mw, vz_mw
