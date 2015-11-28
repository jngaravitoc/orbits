import numpy as np
import sys
from parameters import *
from acceleration import *
from leapfrog import *
from output import *


t, xsat , ysat, zsat, xhost, yhost, zhost, vxsat, vysat, vzsat, vxhost, vyhost, vzhost = leapfrog()

writing(t, xsat, ysat, zsat, vxsat, vysat, vzsat, xhost, yhost, zhost, vxhost, vyhost, vzhost)
