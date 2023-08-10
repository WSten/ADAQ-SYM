# ADAQ-SYM
ADAQ-SYM is a code that performs symmetry analysis of defect orbitals,
finding irreducible representations (IR) of orbitals and symmetry allowed transitions between orbitals.
The code only works for defects simulated in VASP.

## Installation

### Requirements

* numpy
* matplotlib
* AFLOW-SYM, follow the instructions at http://www.aflowlib.org/install-aflow/, add AFLOW to your PATH and PYTHONPATH
* VaspBandUnfolding, follow the instructions at https://github.com/QijingZheng/VaspBandUnfolding
* Character table files, see instructions below http://gernot-katzers-spice-pages.com/character_tables/

### Instructions

To download this repository use:
```bash
git clone https://github.com/WSten/ADAQ-SYM.git
wget http://www.gernot-katzers-spice-pages.com/character_tables/character_tables.tar.xz
mkdir character_tables
tar -xf character_tables.tar.xz -C character_tables
sed -i "s|/path/to|$PWD|g" ADAQ-SYM/settings.json
```
(Make this a python module so that its importable)

### Example
To try ADAQ-SYM on a NV(-) center do the following:
```bash
wget ......dedur01/data/public/ADAQ-SYM/NV_example
cd NV_example
python3 /path/to/ADAQ-SYM/run_main.py 1022 1024 1022 1024
```

## Use

Starting with a completed simulation of a defect in VASP, the ```main()``` function in ```analysis.py``` reads the outputs from VASP and performs the symmetry analysis by first finding characters by calculating overlaps, and then assigning irreducible representations to each orbital. Then, selection rules are applied to the transition dipole moment to find allowed transitions.
The continuous symmetry measure (CSM) is also evaluated to see how well orbital conform to different irreducible representations.

### Required Inputs
The files required from a vasp simulation are:
Wave function in ```WAVECAR```, atomic positions in ```POSCAR``` or ```CONTCAR```, and eigenvalues and occupation in ```EIGENVAL```.
The default names of these files the code looks for are: ```WAVECAR```, ```CONTCAR``` and ```EIGENVAL```. If you have other filenames they need to be specified as arguments to main().

Two arrays with indices of the considered bands in each spin channel are required.
E.g. ```main([1025,1026,1127], [1024,1025,1026,1127])```. This can easily be generated with the ```run_main.py``` script.

There are options of setting a nametag to the output files, and where the output files are written.

### Settings
Settings of the code can be found and edited in the ```settings.json``` file.
The only setting required to be set by the user points to the directory with character tables.
```"char_table_dir": "/path/to/character_tables"```

Other settings can be left as default, and changed when required, see "article".

### Overlap Calculations
For each orbital the center of mass is calculated, and this centers is used a the fixed point when the symmetry operators U are applied to the wave function  of orbital i, ```psi_i```. The overlap, or symmetry operator expectation value (SOEV) is calculated as:
```<psi_i | U | psi_i>```.


### Symmetry Analysis
The character of a class is taken to be the average of the SOEVs of that class.
These characters are then projected on an IR 'Gamma', and if this projection is within the "IR_tolerance" of 1, then the orbital is said to transform as IR 'Gamma'.
This is done for all IRs in the relevant point group.
CSM is calculated based on this projection to give a numerical measure of how well orbitals conforms to different IRs.

### Outputs
The centers are written to files: ```Centers_*_Sx.npy``` and ```Centers_*_Sx.txt```.
The overlaps are written to files: ```Overlaps_*_Sx.pickle``` and ```Overlaps*.txt```.
The characters, IRs and transitions are written to files: ```Transitions_*_Sx.pickle``` and ```Transitions_*_Sx.txt```.
CSM is written to ```CSM*.txt```.

'*' is a name chosen by the user, and 'x' = 1 or 2 for spin up or down, respectively.

### Plots
```plot_levels()``` in ```plotting.py``` generates an eigenvalue level diagram of the single particle states,
where irreducible representation is shown for each level, and allowed transitions
are marked by arrows with color depending on polarization.

## Cite
Cite as ...
https://arxiv.org/abs/2307.04381
