import numpy as np
import sys
from parameters import *
from acceleration import *
from leapfrog import *

t, xx , yy, zz, xx_mw, yy_mw, zz_mw, vx, vy, vz, vx_mw, vy_mw, vz_mw = leapfrog()


#R = np.sqrt((xx-xx_mw)**2 + (yy-yy_mw)**2 + (zz-zz_mw)**2)
R = np.sqrt((xx)**2 + (yy)**2 + (zz)**2)
RMW = np.sqrt(xx_mw**2 + yy_mw**2 + zz_mw**2)

print "# positions and velocities are relative to a (0,0,0) point"
print "# t (Gyr), x_sat(kpc), y_sat(kpc), z_sat(kpc), xmw(kpc), ymw(kpc), zmw(kpc), vx_sat(km/s), vy_sat(km/s), vz_sat(km/s) ,  vxmw(km/s), vymw(km/s), vzmw(km/s)"

for i in range(len(xx)):
    print t[i], xx[i], yy[i], zz[i], xx_mw[i], yy_mw[i], zz_mw[i],
vx[i] / conv_factor,  vy[i] / conv_factor , vz[i] / conv_factor ,
vx_mw[i] / conv_factor, vy_mw[i] / conv_factor, vz_mw[i] / conv_factor

