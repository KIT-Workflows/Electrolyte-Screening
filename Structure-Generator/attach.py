import ase
from ase.io import write, read
from ase import Atoms
from ase.build import attach 
import numpy as np
import os, sys, tarfile, yaml


########### Redaing inputs ############
with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)

tar_file = wano_file["structures"]
tar = tarfile.open(tar_file)
tar.extractall()
tar.close()


var_len = len(wano_file["Add-Molecules"])
mol_names = []
distances = []

for ii in range(var_len):
    mol_names.append('structures/' + wano_file["Add-Molecules"][ii]["mol-name"])
    distances.append(wano_file["Add-Molecules"][ii]["distance"])

print(mol_names)
print(distances)

# mol1_file = read("structures/LiFSA.xyz")
# #write('LiFSA.pdb', mol_file)
# mol2_file = read("structures/CH3CN.xyz")
#write('CH3CN.pdb', mol_file)

for ii in range(var_len-1):
    if ii == 0:
        mol_att = attach.attach_randomly(read(mol_names[ii]), read(mol_names[ii+1]), distances[ii], rng=np.random)
        write('add_mol.xyz', mol_att)
    else:
        mol2_file = read("add_mol.xyz")
        mol_att = attach.attach_randomly(mol2_file, read(mol_names[ii]), distances[ii], rng=np.random)
        write('add_mol.xyz', mol_att)
        #mol4 = attach.attach_randomly(mol2_file,mol3, 2.0, rng=np.random)
    
