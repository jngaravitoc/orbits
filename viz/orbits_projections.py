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

#t_sat = data_sat[:,0]
x_sat = data_sat[:,1]
y_sat = data_sat[:,2]
z_sat = data_sat[:,3]
#vx_sat = data_sat[:,4]
#vy_sat = data_sat[:,5]
#vz_sat = data_sat[:,6]

#t_host = data_host[:,0]
x_host = data_host[:,1]
y_host = data_host[:,2]
z_host = data_host[:,3]
#vx_host = data_host[:,4]
#vy_host = data_host[:,5]
#vz_host = data_host[:,6]


'''
plt.plot(y_sat, z_sat, lw=2)
plt.scatter(y_sat[0], z_sat[0], s=60)
plt.plot(y_host, z_host, ls='--', lw=2, c='k')
plt.scatter(y_host[0], z_host[0], s=60, c='k')
plt.ylabel("y(kpc)", fontsize=25)
plt.xlabel("x(kpc)", fontsize=25)
plt.show()
'''
plt.plot(x_sat, z_sat, lw=2)
plt.plot(x_host, z_host, ls='--', lw=2, c='k')
plt.scatter(x_sat[0], z_sat[0], s=60)
plt.scatter(x_host[0], z_host[0], s=60, c='k')
plt.ylabel("z(kpc)", fontsize=25)
plt.xlabel("x(kpc)", fontsize=25)
plt.show()

