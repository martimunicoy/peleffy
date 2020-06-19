# Global imports
from collections import defaultdict


class Rotamer(object):
    def __init__(self, atom1, atom2, resolution=30):
        self.atom1 = atom1
        self.atom2 = atom2
        self.resolution = resolution

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def resolution(self):
        return self._resolution


class RotamerLibrary(object):
    def __init__(self, residue_name):
        self._residue_name = residue_name
        self._rotamers = defaultdict(list)

    def add_rotamer(self, rotamer, group_id):
        self._rotamers[group_id].append(rotamer)

    def to_file(self, path):
        with open(path, 'w') as file:
            file.write('rot assign res {} &\n')
            for i, group in enumerate(self.rotamers.keys()):
                if i > 0:
                    file.write('     newgrp &\n')
                for rotamer in self.rotamers[group]:
                    file.write('   sidelib FREE{} {} {} &\n'.format(
                        rotamer.atom1, rotamer.atom2, rotamer.resolution))

    @property
    def residue_name(self):
        return self._residue_name

    @property
    def rotamers(self):
        return self._rotamers


class RotamerLibraryBuilder(object):
    def build_from_molecule(self, molecule):
        try:
            from rdkit import Chem
        except ImportError:
            raise Exception('RDKit Python API not found')
