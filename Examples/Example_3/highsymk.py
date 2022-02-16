from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.inputs import Incar

import shutil
import os
import numpy as np

from pymatgen.core.structure import IStructure

class vaspinput(object):
    '''
    generate INCAR parameter
    '''
    def __init__(self, structure, todo='self', soc=False,mag=False,saxis=False,nband=True):
        self.todo = todo
#        structure = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
        self.structure = SpacegroupAnalyzer( structure ).get_primitive_standard_structure()
        self.structure = structure
#        structure = SpacegroupAnalyzer(structure).find_primitive()
        Poscar(self.structure).write_file('POSCAR')
	

    def kpoints(self, denband=30):
        todo = self.todo
        if todo == 'band':
            path = HighSymmKpath(self.structure)
            kfile = Kpoints.automatic_linemode(denband, path)
            kfile.write_file('KPOINTS')


if __name__ == '__main__':
    structure=IStructure.from_file('POSCAR')
    a = vaspinput(structure, todo='band')
    a.kpoints()
