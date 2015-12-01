import numpy as np
import matplotlib.pyplot as plt

sat = [3, 5, 8, 10, 18, 25]

plt.figure(figsize=(6, 9))
plt.subplot(2, 1, 1)
for i in sat:
     data_host = np.loadtxt('../data/LMC_'+str(i)+'E10_1E12_free_host.txt')
     data_sat = np.loadtxt('../data/LMC_'+str(i)+'E10_1E12_free_sat.txt')
     Rgal = np.sqrt((data_sat[:,1]-data_host[:,1])**2 + (data_sat[:,2] -
     data_host[:,2])**2 + (data_sat[:,3] - data_sat[:,3])**2)
     time = data_host[:,0]
     plt.plot(time, Rgal, lw=2)
plt.xlabel('$\mathrm{time(Gyr)}$', fontsize=18)
plt.ylabel('$\mathrm{R_{gal}(Kpc)}$', fontsize=18)

plt.subplot(2, 1, 2)

for i in sat:
     data_host2 = np.loadtxt('../data/LMC_'+str(i)+'E10_2E12_free_host.txt')
     data_sat2 = np.loadtxt('../data/LMC_'+str(i)+'E10_2E12_free_sat.txt')
     Rgal = np.sqrt((data_sat2[:,1]-data_host2[:,1])**2 + (data_sat2[:,2] -
     data_host2[:,2])**2 + (data_sat2[:,3] - data_sat2[:,3])**2)
     time2 = data_host2[:,0]
     plt.plot(time2, Rgal, lw=2)
plt.xlabel('$\mathrm{time(Gyr)}$', fontsize=18)
plt.ylabel('$\mathrm{R_{gal}(Kpc)}$', fontsize=18)
plt.show()
