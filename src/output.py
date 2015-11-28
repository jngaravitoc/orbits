from parameters import *

def writing(t, xsat, ysat, zsat, vxsat, vysat, vzsat, xhost, yhost, zhost, vxhost, vyhost, vzhost):

    satellite = open(path + filename + "_sat.txt", 'w')
    host = open(path + filename + "_host.txt", 'w')
    print "writing orbits in:", path
    for i in range(len(t)):
        satellite.write("%f %f %f %f %f %f %f \n"%(t[i], xsat[i], ysat[i], zsat[i], vxsat[i], vysat[i], vzsat[i]))
        host.write("%f %f %f %f %f %f %f \n"%(t[i], xhost[i], yhost[i], zhost[i], vxhost[i], vyhost[i], vzhost[i]))
    satellite.close
    host.close
    print "done!"
