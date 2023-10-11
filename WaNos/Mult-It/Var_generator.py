import numpy as np
import yaml, sys, glob
import shutil, os, tarfile
from os import listdir
from pathlib import Path
from os.path import isfile, join

def sub_directory(folder):
    mylist=[]
    for root,subdirectories, files in os.walk(folder):
        for subdirectory in subdirectories:
            sub=os.path.join(root, subdirectory)
    for file in files:
        list_str=os.path.join(root, file)
        new_string = list_str.replace(folder, "")
        mylist.append(new_string[1:])
    return mylist


def get_filepaths(directory):
    
    '''
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    '''

    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


def extract_tar(tar_file):
    path = os.path.abspath(os.getcwd())
    opened_tar = tarfile.open(tar_file)
     
    if tarfile.is_tarfile(tar_file):
        opened_tar.extractall(path)
    else:
        print("The tar file you entered is not a tar file")

def tardir(path, tar_name):
       with tarfile.open(tar_name, "w:gz") as tar_handle:
              for root, dirs, files in os.walk(path):
                     for file in files:
                            tar_handle.add(os.path.join(root, file))


if __name__ == '__main__':
       
       float_outdict = {}
       float_outdict["iter"] = []
       with open("float_output_dict.yml",'w') as out:
              yaml.dump(float_outdict, out,default_flow_style=False) 

       int_outdict = {}
       int_outdict["iter"] = []
       with open("int_output_dict.yml",'w') as out:
              yaml.dump(int_outdict, out, default_flow_style=False)

       structure_outdict = {}
       structure_outdict["iter"] = []

       with open("structure_output_dict.yml",'w') as out:
              yaml.dump(structure_outdict, out,default_flow_style=False)

       with open("output_dict.yml",'w') as out:
              yaml.dump(structure_outdict, out, default_flow_style=False)

       
       ''' input paprameters'''
       
       with open("rendered_wano.yml") as file:
              wano_file = yaml.full_load(file)
       
       float_var = wano_file['Float']
       float_int = wano_file['Int']
       file_struct = wano_file['Structures']

       ''' Writing the float output file '''
       
       if float_var:
              varF_begin = wano_file['VarF-begin']
              varF_end = wano_file['VarF-end']
              n_points = wano_file['N-points']
              float_outdict = {}
              float_outdict["iter"] = []
              temp_var = np.linspace(varF_begin, varF_end, n_points)
              float_outdict["iter"] = [float(i) for i in temp_var]
              print(float_outdict)

              with open("float_output_dict.yml",'w') as out:
                     yaml.dump(float_outdict, out,default_flow_style=False)

              shutil.copy2("float_output_dict.yml", "output_dict.yml")

       if float_int:
              varI_begin = wano_file['VarF-begin']
              varI_end = wano_file['VarF-end']
              step_var = wano_file['Step']
              int_outdict = {}
              int_outdict["iter"] = []
              temp_var = np.arange(varI_begin, varI_end, step_var)
              int_outdict["iter"] = [int(i) for i in temp_var]

              with open("int_output_dict.yml",'w') as out:
                     yaml.dump(int_outdict, out, default_flow_style=False)

              shutil.copy2("int_output_dict.yml", "output_dict.yml")

       if file_struct:
              if os.path.isdir('Structures'):
                     shutil.rmtree('Structures')
              
              if os.path.isfile("Structures.tar"):
                     fname = "Structures.tar"
                            
                     if fname.endswith("tar.gz"):
                            tar = tarfile.open(fname)
                            tar.extractall("Structures")
                            tar.close()
                     elif fname.endswith("tar.xz"):
                            tar = tarfile.open(fname)
                            tar.extractall("Structures")
                            tar.close()    
                     elif fname.endswith("tar"):
                            tar = tarfile.open(fname)
                            tar.extractall("Structures")
                            tar.close()
              folder=os.path.join(os.getcwd(),"Structures")

              
              my_list = sub_directory(folder)
              structure_outdict = {}
              structure_outdict["iter"] = my_list

              with open("structure_output_dict.yml",'w') as out:
                     yaml.dump(structure_outdict, out,default_flow_style=False)

              structure_outdict["struct_len"] = len(structure_outdict["iter"])
              with open("output_dict.yml",'w') as out:
                     yaml.dump(structure_outdict, out, default_flow_style=False)

       else:
              print ("Structure dir does not exist. Creating an empty dir")
              tar = tarfile.open("Structures.tar", "w")
              os.mkdir('Structures') 
              tardir('./Structures', 'Structures.tar')
              tar.close()
              shutil.rmtree('Structures')
              var_str = []

              structure_outdict = {}
              structure_outdict["iter"] = var_str               