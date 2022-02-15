# Automatic-wannier-flow
Automatic construction of wannier functions for any 3D transition metal based system with or without SOC

Written by Ilias Samathrakis and Zeying Zhang

*****************************************************************
INTRODUCTION

We present an open-source code, written in Python, that is used to
construct wannier functions for any 3D transition metal based system
with or without SOC automatically. 

*****************************************************************
REQUIREMENTS

The code is only tested on Linux systems and requires Python3 and the library 'pymatgen' 

Download from:
https://www.python.org/downloads/
and
https://pymatgen.org/

The code requires the following external software: 
1) VASP -> tested in version 5.2.12
2) wannier90 -> tested in version 2.0.0

Download from:
http://www.wannier.org/download/

*****************************************************************
DESCRIPTION OF FILES

1) vasp_input_hte_cif.py -> source code
2) input_parameters.py   -> input of source code
3) autoconstruction.py   -> constructs wannier90 input file
4) compareband.py        -> compares band structures obtained from DFT and WFs
5) highsymk.py           -> generates high symmetric kpaths for the KPOINTS file 
6) plotband.py           -> plots the band structure obtained from VASP
7) Notes.pdf             -> Detailed explanation

*****************************************************************
HOW TO USE

1) Paste the files within the directory you want to perform the calculation.
2) POSCAR file suitable for VASP calculation should be included in the directory.
3) Modify the file 'input_parameters.py' based on your system
4) Modify the lines 344-362 of the file 'vasp_input_hte_cif.py' based on your operating system

For detailed explanation read Notes.pdf

*****************************************************************
OUTPUT

Run the code using the command
python vasp_input_hte_cif.py

After executing the code, the files needed for the complete calculation and the submission script are generated.

Inspect and submit!
