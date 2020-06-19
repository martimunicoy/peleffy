
# Global imports
import offPELE


class Impact(object):
    def __init__(self, molecule):
        if isinstance(molecule, offPELE.topology.Molecule):
            self._initialize_from_molecule(molecule)
        else:
            raise Exception('Invalid input molecule for Impact template')

    def _initialize_from_molecule(self, molecule):
        self._molecule = molecule

    def write(self, path):
        with open(path) as file:
            self._write_header(file)
            self._write_resx(file)
            self._write_nbon(file)
            self._write_bond(file)
            self._write_thet(file)
            self._write_phi(file)
            self._write_iphi(file)
            self._write_end(file)

    def _write_header(self, file):
        file.write('* LIGAND DATABASE FILE ({})\n'.format(
            self.molecule.forcefield))
        file.write('* File generated with offPELE {}\n'.format(
            offPELE.__version__))

    def _write_resx(self, file):
        pass

    def _write_nbon(self, file):
        pass

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
