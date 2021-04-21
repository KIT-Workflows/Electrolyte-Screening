import numpy as np
import os, sys, re, yaml

def check_vdw(dict_wano):
    temp_var = ""
    vdw_var = ""
    if dict_wano["Approach"]["Dispersion-corrections"]:
        temp_var = 'dsp\n'
        if dict_wano["Approach"]["vdw-corrections"] == "disp2":
            vdw_var = "old\n"
        elif dict_wano["Approach"]["vdw-corrections"] == "disp3":
            vdw_var = "on\n"
        elif dict_wano["Approach"]["vdw-corrections"] == "disp3 -bj":
            vdw_var = "bj\n"
        else:
            temp_var = "\n"
            vdw_var = "\n"
    
    return  temp_var, vdw_var


def check_optmizer(dict_wano):
    optmizer_var = str()
    if dict_wano["Approach"]["Optimizer"] == "single-shot":
        optmizer_var = optmizer_var
    else:
        optmizer_var = "jobex -dscf -c 3000 > jobex.out\n"

    return optmizer_var

if __name__ == '__main__':
         
    with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)

    #print(wano_file["Molecular-structure"])
    # dict_define = get_turbomole_inputs()

    # # Reading the inputs from param file
    # for var_key, var_value in wano_file["Approach"].items():
    #     if var_key in dict_define.keys():
    #         dict_define[var_key] = var_value
    #         print(var_key, var_value)
    #     else:
    #         dict_define[var_key] = var_value

    # print(len(dict_define))

    file_name = "define_file.sh"
    with open(file_name, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('\n')
        
        f.write('module purge\n')
        f.write('module load turbomole/7.4.1\n')
        f.write('\n')

        f.write('x2t ' + wano_file["Molecular-structure"]["Structure-file"] + ' > coord'+ '\n')
        f.write('\n')

        f.write('define << EOF > define.out\n')
        f.write('\n')
        f.write('\n')
        f.write('\n')

        f.write('a coord\n')
        f.write('*\n')
        f.write('\n')

        if wano_file["Molecular-structure"]["Symmetry"]:
            f.write('desy'+ ' ' + str(wano_file["Molecular-structure"]["Threshold"]) + '\n')
            f.write('\n')

        if wano_file["Molecular-structure"]["Internal-coordinates"]:
            f.write('ired\n')
            f.write('*\n')
            f.write('\n')
        
        else:
            f.write('*\n')
            f.write('no\n')
            f.write('\n')

        f.write('b all'+ ' ' + str(wano_file["Basis-set"]["Basis-set-types"]) + '\n')
        f.write('*\n')
        f.write('eht\n')
        f.write('\n')
        f.write(str(wano_file["Starting-orbitals"]["Charge"])+'\n') 
        f.write('\n')
        f.write('\n')
        # # new added
        # f.write('scf\n')
        # f.write('iter\n')
        # f.write('2000\n')
        # f.write('damp\n')
        # f.write('20\n')        
        # f.write('\n')
        # f.write('\n')
        #############
        f.write('ri\n')
        f.write('on\n')
        f.write('m 5000\n')
        f.write('*\n')
        # Method
        f.write(str(wano_file["Approach"]["Method"])+'\n')
        f.write('on\n')
        f.write('func ' + str(wano_file["Approach"]["Functional"]) + '\n')
        f.write('grid m4\n')
        f.write('*\n')
        # check vdw correction
        print(check_vdw(wano_file)[0])
        print(check_vdw(wano_file)[1])

        f.write(check_vdw(wano_file)[0])
        f.write(check_vdw(wano_file)[1])
        f.write('\n')
        
        f.write('*\n')

        f.write('\n')
        f.write('\n')
        
        f.write('EOF\n')
        f.write('\n')

        f.write('ridft > ridft.out\n')
        f.write('eiger > eiger.out\n')
        f.write('\n')

        f.write(check_optmizer(wano_file))
        f.write('\n')

        f.write('exit 0\n')

    os.system("chmod +x " + file_name)
    print("define file successfully created")