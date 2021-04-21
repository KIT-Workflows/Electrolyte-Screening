#!/home/ws/gt5111/testEnv/bin/python
import numpy as np
import sys, os

try:
  import yaml
except ImportError:
  print ("Trying to Install required module: requests\n")
  os.system('pip3 install PyYAML')
import yaml



## inputs #############################
with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)
#######################################

outdict = {}
outdict["iter"] = []
if wano_file["Float"]:
    var_begin = wano_file["VarF-begin"]
    var_end = wano_file["VarF-end"]
    var_npoints = wano_file["N-points"]
    
    step = (var_end-var_begin)/var_npoints
    for i in range(0,var_npoints):
        z_0 = var_begin + step*i
        z_0 = round(z_0,3)
        outdict["iter"].append(float(z_0))    
elif wano_file["Int"]:
    var_begin = wano_file["VarI-begin"]
    var_end = wano_file["VarI-end"]
    step = wano_file["Step"]
    if step == 0:
        print( "N_points must fit in the interval" )
        sys.exit()
        
    for z_0 in range(var_begin,var_end, step):
        #z_0 = var_begin + int(step*i)
        outdict["iter"].append(z_0)

with open("output_dict.yml",'w') as out:
    yaml.dump(outdict, out,default_flow_style=False)

print(type(f_file["iter"]))
# with open("output_dict.yml",'w') as out:
#     yaml.safe_dump(outdict, out)