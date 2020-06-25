
# Global imports
import offPELE


class Impact(object):
    def __init__(self, molecule):
        if (isinstance(molecule, offPELE.topology.Molecule)
                or isinstance(molecule, offPELE.topology.molecule.Molecule)):
            self._initialize_from_molecule(molecule)
        else:
            raise Exception('Invalid input molecule for Impact template')

    def _initialize_from_molecule(self, molecule):
        self._molecule = molecule

    def write(self, path):
        with open(path, 'w') as file:
            self._write_header(file)
            self._write_resx(file)
            self._write_nbon(file)
            self._write_bond(file)
            self._write_thet(file)
            self._write_phi(file)
            self._write_iphi(file)
            self._write_end(file)

    def _write_header(self, file):
        file.write('* LIGAND DATABASE FILE')
        if self.molecule.forcefield:
            file.write(' ({})'.format(self.molecule.forcefield))
        file.write('\n')
        file.write('* File generated with offPELE {}\n'.format(
            offPELE.__version__))

    def _write_resx(self, file):
        # template name
        file.write('{:5}'.format(self.molecule.name))
        # number of non bonding parameters
        file.write('{:6d}'.format(len(self.molecule.atoms)))
        # number of bond parameters
        file.write('{:6d}'.format(len(self.molecule.bonds)))
        # number of angle parameters
        file.write('{:7d}'.format(len(self.molecule.angles)))
        # number of dihedral parameters
        # TODO doublecheck that it is indeed the sum of propers and impropers
        file.write('{:7d}'.format(len(self.molecule.propers)
                                  + len(self.molecule.impropers)))
        # # number of non-null elements in the interactions matrix
        # TODO It might not be always 0
        file.write('{:8d}'.format(0))
        file.write('\n')

    def _write_nbon(self, file):
        for atom in self.molecule.atoms:
            # atom id number
            file.write('{:5d}'.format(atom.index))

    def _write_bond(self, file):
        pass

    def _write_thet(self, file):
        pass

    def _write_phi(self, file):
        pass

    def _write_iphi(self, file):
        pass

    def _write_end(self, file):
        file.write('END\n')

    @property
    def molecule(self):
        return self._molecule
