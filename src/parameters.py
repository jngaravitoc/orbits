#########################################################################
# This code reads the input parameters the code uses the function and
# method described in:
# https://wiki.python.org/moin/ConfigParserExamples
#########################################################################
import ConfigParser
import sys

input_param = sys.argv[1]

Config = ConfigParser.ConfigParser()
Config.read(input_param)


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

x_sat = float(ConfigSectionMap("params")['xsat'])
y_sat = float(ConfigSectionMap("params")['ysat'])
z_sat = float(ConfigSectionMap("params")['vzsat'])
vx_sat = float(ConfigSectionMap("params")['vxsat'])
vy_sat = float(ConfigSectionMap("params")['vysat'])
vz_sat = float(ConfigSectionMap("params")['vzsat'])
M_sat = float(ConfigSectionMap("params")['msat']) * 1E10
Rvir_sat = float(ConfigSectionMap("params")['rvirsat'])
c_sat = float(ConfigSectionMap("params")['csat'])
x_host = float(ConfigSectionMap("params")['xhost'])
y_host = float(ConfigSectionMap("params")['yhost'])
z_host = float(ConfigSectionMap("params")['zhost'])
vx_host = float(ConfigSectionMap("params")['vxhost'])
vy_host = float(ConfigSectionMap("params")['vyhost'])
vz_host = float(ConfigSectionMap("params")['vzhost'])
M_host = float(ConfigSectionMap("params")['mhost']) * 1E10
Rvirhost =  float(ConfigSectionMap("params")['rvirhost'])
chost = float(ConfigSectionMap("params")['chost'])
M_disk = float(ConfigSectionMap("params")['mdisk']) * 1E10
a_disk = float(ConfigSectionMap("params")['adisk'])
b_disk = float(ConfigSectionMap("params")['bdisk'])
M_bulge = float(ConfigSectionMap("params")['mbulge']) * 1E10
rh_disk = float(ConfigSectionMap("params")['rhdisk'])
Host_move = float(ConfigSectionMap("params")['hostmove'])
alpha_df = float(ConfigSectionMap("params")['alpha'])
Host_df = float(ConfigSectionMap("params")['hostdf'])
