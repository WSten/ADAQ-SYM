# ADAQ-SYM
What does this code do?

## Installation 

### This Repository

### AFLOW-SYM
Follow the instructions at http://www.aflowlib.org/install-aflow/

### VaspBandUnfolding
Follow the instructions at https://github.com/QijingZheng/VaspBandUnfolding

### Character tables
Download from http://gernot-katzers-spice-pages.com/character_tables/

## Use

main() 

### Required Inputs
The files required from a vasp simulation are: 
Wave function in WAVECAR, atomic positions in POSCAR or CONTCAR, and eigenvalues and occupation in EIGENVAL.
The default names of these files the code looks for are: WAVECAR, CONTCAR and EIGENVAL. If you have other filenames they need to be specified as arguments to main().

Two arrays with indices of the considered bands in each spin channel are required.
E.g. main([1025,1026,1127], [1024,1025,1026,1127]). This can easily be generated with the run_main script.

There are options of setting a nametag to the output files, and where the output files are written.

### Settings
Settings of the code can be found and edited in the "settings.json" file.
The only setting required to be set by the user is point to the directory with character tables.
"char_table_dir": "/path/to/character_tables"

Other settings can be left as default, and changed when requred, see "article".

### Overlap Calculations
For each orbital the center of mass is calculated, and this centers is used a the fixed point when the symmetry operators U are applied to the wave function  of orbital i, psi_i. The overlap, or symmetry operator expecation value (SOEV) is calculated as:
<psi_i | U | psi_i>

The centers are written to files: "Centers_*_Sx.npy" and "Centers_*_Sx.txt".
The overlaps are written to files: "Overlaps_*_Sx.pickle" and "Overlaps_*_Sx.pickle"

### Symmetry Analysis


### Outputs 


### Plots


## Licence
Some open-source licence 

## Cite
Cite as ...



