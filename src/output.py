from parameters import *

def writing(t, xsat, ysat, zsat, vxsat, vysat, vzsat, xhost, yhost, zhost, vxhost, vyhost, vzhost):
    satellite = open(path + filename + "_sat.txt", 'w')
    host = open(path + filename + "_host.txt", 'w')
    print len(xhost), len(t)
    print "writing orbits in:", path
    for i in range(len(t)):
        satellite.write("%f %f %f %f %f %f %f \n"%(t[i], xsat[i],
ysat[i], zsat[i], vxsat[i]/conv_factor, vysat[i]/conv_factor,
vzsat[i]/conv_factor))
        host.write("%f %12f %12f %12f %12f %12f %12f \n"%(t[i],
xhost[i], yhost[i], zhost[i], vxhost[i]/conv_factor,
vyhost[i]/conv_factor, vzhost[i]/conv_factor))
    satellite.close
    host.close
    print "done!"
