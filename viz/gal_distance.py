import numpy as np
import matplotlib.pyplot as plt
import sys

if (len(sys.argv)<2):
    print "use: satellite  host"
filename1 = str(sys.argv[1])
filename2 = str(sys.argv[2])
path1 = "../data/" + filename1
path2 = "../data/" + filename2

data_sat = np.loadtxt(path1)
data_host = np.loadtxt(path2)

t_sat = data_sat[:,0]
x_sat = data_sat[:,1]
y_sat = data_sat[:,2]
z_sat = data_sat[:,3]
vx_sat = data_sat[:,4]
vy_sat = data_sat[:,5]
vz_sat = data_sat[:,6]

t_host = data_host[:,0]
x_host = data_host[:,1]
y_host = data_host[:,2]
z_host = data_host[:,3]
vx_host = data_host[:,4]
vy_host = data_host[:,5]
vz_host = data_host[:,6]

R_sat = np.sqrt(x_sat**2 + y_sat**2 + z_sat**2)
R_host = np.sqrt(x_host**2 + y_host**2 + z_host**2)
R_gal = np.sqrt((x_sat - x_host)**2 + (y_sat - y_host)**2 + (z_sat - z_host)**2)

plt.plot(t_sat, R_gal, lw=2)
#plt.plot(t_sat, R_sat, ls='--', lw=2, c='k')
plt.ylabel("Galactocentric Radius", fontsize=25)
plt.xlabel("time(Gyrs)", fontsize=25)
plt.show()
