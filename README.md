# Electrolyte-Screening

In this workflow, we use the SimStack framework features to screening electrolytes systems using DFT calculation. Here, we combine four different WaNos: Range-It, Structure-Generator, DFT-Turbomole, and Table-Generator, to set up an electrolyte system, load the file structure, and choose the methods embedded in the DFT approach using Turbomole code. A table containing the system's HOMO-LUMO gap energy and molecule label is the expected output of this protocol.

Using the drag-and-drop environment of Simstack, we can build the workflow depicted in **Fig 1** in four steps. The Range-It WaNo accounts for the number of different configurations of a given system.  In the second step, we add the Structure-Generator **WaNo** inside the ForEach loop control to generate the ```.xyz``` files of the configuration system. In the third step, we insert the DFT-Turbomole **WaNo**, which will receive the generated files from the previous one. At this step, We can take advantage of the parallelization in the HPC remote resources once the ForEach loop control is designed for this end.  Table-Generator **WaNo** extracts two variable values on the ```job.last``` file: the output file of steps two and three. This **WaNo** builds a table named Table_var in CSV format at the end of the protocol.

### In this workflow, we will be able to:
```
1. To set up many electrolyte configurations from an initial seed (Range-It).
2. Load a molecule seed and attach many other molecules to the seed (Structure-Generator).
3. Run the geometric DFT calculations using Turbomole code, accounting for the proper corrections (DFT-Turbomole).
4. Arrange all the HOMO-LUMO gap energy values of the system in a table format (Table-Generator).
```

## DFT_Surface workflow with **_ForEach_** loop control
![Semantic description of image](Electrolyte-Screening.png)

**Fig 1** _This workflow aims to perform several DFT calculations of electrolyte systems. It comprises Range-It, Structure-Generator, DFT-Turbomole, and Table-Generator WaNos connected by the ForEach loop control. In step 1, we generate the number of configurations. Steps 2 and 3 define the electrolyte designs and the DFT calculation methods employed in the simulation. The **WaNo** in the last step extracts the inquired variables of the output file from the previous actions._

aqui

## 1. Python Setup
To get this workflow up running on your available computational resources, make sure to have the below libraries installed on Python 3.6 or newer.

```
1. Atomic Simulation Environment (ASE).
2. Python Materials Genomics (Pymatgen).
3. Numpy, os, sys, re, yaml, subprocess.
4. json, csv, shutil, tarfile. 
```
## 2. Mult_Mol Inputs
- Range of the variable position. 
- Number of points in the present in the range. 
- Beginning of the Molecule name, which should appear in all molecules. 
- Directory with the zip file of the molecules.
## 3. Mult_Mol Output
- It should pass all the information to the next WaNo inside the ForEach loop through the `Mult_Mol.iter.*` command on the top of the loop, as Fig 1 shows in step 2.
## 4. Surface Inputs
- Aux_var should be set as `${ForEach_iterator_ITER}` from import workflow variable.
- Mol_name should be set as `Mult_Mol.Molecule_name` from import workflow variable.
- Defining bulk unit cell types, element, and lattice constant.
- Defining slab size, vacuum size, Miller index of the surface, and as an option, set a supercell.
- Check the box when you want to adsorb a molecule on the surface previously defined.
- Setting the molecule distance over the surface and molecule-molecule image distance.
## 5. Surface Output
- POSCAR and Input_data.yml files, which should be passed to DFT_VASP WaNo.
## 6. DFT_VASP Inputs
- **INCAR tab**: as an option, we can set all INCAR flags available within VASP. However, we expose only a few of them, which are essential for the problem. See the GUI of this WaNo. A brief description of each flag pops up when we rover the mouse over the inputs.
- **KPOINTS tab**: Here the user can define two types of KPOINTS, `Kpoints_length` and `Kpoints_Monkhorst`.
- **Analysis tab**: Aimed to compute Bader charge analysis and DOS.
- **Files_Run tab**: Mandatory loads the POSCAR file, and as an option can load INCAR, POTCAR, KPOINTS, and KORINGA files. The KORINGA file can be any file. In the case of this problem, it loads the Input_data.yml file.
## 7. DFT_VASP Output
- OUTCAR file.
## 8. Table_Generator Inputs
- Search_in_File: Should be set as OUTCAR and import the OUTCAR file using `ForEach/*/DFT_VASP/outputs/OUTCAR` command.
- Delete_Files: check the box option.
- Search_Parameters: Set the variables `z_0`, `File_number`, and `energy`.  
## 8. Table_Generator Output
- Table_var file in CSV format containing the variables defined in the Search_Parameters field.
