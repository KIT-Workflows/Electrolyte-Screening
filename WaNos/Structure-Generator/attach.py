import ase
from ase.io import write, read
from ase import Atoms
from ase.build import attach 
import numpy as np
import os, sys, tarfile, yaml

def load_yaml_file(filepath):
    """Load a YAML file and return its content."""
    with open(filepath) as file:
        return yaml.full_load(file)

def extract_tar_file(tar_filepath):
    """Extract a tar file to the current directory."""

    folder_names = []
    with tarfile.open(tar_filepath, 'r') as tar:
        folder_names.extend(
            member.name for member in tar.getmembers() if member.isdir()
        )

    with tarfile.open(tar_filepath) as tar:
        tar.extractall()

    return folder_names

def get_molecule_data_from_yaml(yaml_content, folder_name):
    """Extract molecule names and distances from YAML content."""
    var_len = len(yaml_content["Add-Molecules"])
    mol_names = [f"{folder_name}/" + yaml_content["Add-Molecules"][ii]["mol-name"] for ii in range(var_len)]
    distances = [yaml_content["Add-Molecules"][ii]["distance"] for ii in range(var_len)]
    return mol_names, distances

def attach_molecules(mol_names, distances):
    """Attach molecules randomly and save the result."""
    for ii in range(len(mol_names)-1):
        if ii == 0:
            mol_att = attach.attach_randomly(read(mol_names[ii]), read(mol_names[ii+1]), distances[ii], rng=np.random)
            write('add_mol.xyz', mol_att)
        else:
            mol2_file = read("add_mol.xyz")
            mol_att = attach.attach_randomly(mol2_file, read(mol_names[ii]), distances[ii], rng=np.random)
            write('add_mol.xyz', mol_att)

def main():
    # Load YAML file
    yaml_content = load_yaml_file('rendered_wano.yml')
    
    # Extract tar file
    tar_filepath = yaml_content["structures"]
    folder_names = extract_tar_file(tar_filepath)
    
    # Get molecule names and distances
    mol_names, distances = get_molecule_data_from_yaml(yaml_content,folder_names[0])
    distances = [float(s) for s in distances]
    # print(mol_names)
    print(distances)
    
    # Attach molecules
    attach_molecules(mol_names, distances)

if __name__ == '__main__':
    main()