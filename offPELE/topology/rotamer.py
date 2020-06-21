# Global imports
from collections import defaultdict


class Rotamer(object):
    def __init__(self, atom1, atom2, resolution=30):
        self._atom1 = atom1
        self._atom2 = atom2
        self._resolution = resolution

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
    def __init__(self, residue_name='LIG'):
        self._residue_name = residue_name
        self._rotamers = defaultdict(list)

    def add_rotamer(self, rotamer, group_id):
        self._rotamers[group_id].append(rotamer)

    def to_file(self, path):
        with open(path, 'w') as file:
            file.write('rot assign res {} &\n'.format(self.residue_name))
            for i, group in enumerate(self.rotamers.keys()):
                if i > 0:
                    file.write('     newgrp &\n')
                for rotamer in self.rotamers[group]:
                    file.write('   sidelib FREE{} {} {} &\n'.format(
                        rotamer.resolution, rotamer.atom1, rotamer.atom2))

    @property
    def residue_name(self):
        return self._residue_name

    @property
    def rotamers(self):
        return self._rotamers
