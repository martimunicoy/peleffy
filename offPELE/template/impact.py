
# Global imports
from copy import deepcopy

import offPELE
from offPELE.topology import ZMatrix


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
        file.write('{:6d}'.format(len(self.molecule.angles)))
        # number of dihedral parameters
        # TODO doublecheck that it is indeed the sum of propers and impropers
        file.write('{:8d}'.format(len(self.molecule.propers)
                                  + len(self.molecule.impropers)))
        # # number of non-null elements in the interactions matrix
        # TODO It might not be always 0
        file.write('{:8d}'.format(0))
        file.write('\n')

    def _write_nbon(self, file):
        zmatrix = ZMatrix(self.molecule)

        for i, atom in enumerate(self.molecule.atoms):
            w_atom = WritableAtom(atom)
            # atom id number
            file.write('{:5d}'.format(w_atom.index))
            file.write(' ')
            file.write('{:5d}'.format(w_atom.parent.index))
            file.write(' ')
            file.write('{:1}'.format(w_atom.core))
            file.write('   ')
            file.write('{:4}'.format(w_atom.OPLS_type))
            file.write(' ')
            file.write('{:4}'.format(w_atom.PDB_name))
            file.write(' ')
            file.write('{:5}'.format(w_atom.unknown))
            file.write(' ')
            file.write('{: 11.6f}'.format(zmatrix[i][0]))
            file.write(' ')
            file.write('{: 11.6f}'.format(zmatrix[i][1]))
            file.write(' ')
            file.write('{: 11.6f}'.format(zmatrix[i][2]))
            file.write(' ')
            file.write('\n')

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


class WritableAtom(offPELE.topology.molecule.Atom):
    def __init__(self, atom):
        # We do not want to modify the original object
        atom = deepcopy(atom)

        assert isinstance(atom, (offPELE.topology.molecule.Atom,
                                 offPELE.topology.molecule.DummyAtom)), \
            'Wrong type: {}'.format(type(atom))

        super().__init__(atom.index, atom.core, atom.OPLS_type, atom.PDB_name,
                         atom.unknown, atom.x, atom.y, atom.z, atom.sigma,
                         atom.epsilon, atom.charge, atom.born_radius,
                         atom.SASA_radius, atom.nonpolar_gamma,
                         atom.nonpolar_alpha, atom.parent)

    def none_to_zero(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            if out is None:
                out = int(0)
            return out
        return function_wrapper

    def dummy_to_writable(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            out = WritableAtom(out)
            return out
        return function_wrapper

    def none_to_dummy(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            if out is None:
                out = offPELE.topology.molecule.DummyAtom(index=-1)
            return out
        return function_wrapper

    @property
    @dummy_to_writable
    @none_to_dummy
    def parent(self):
        return super().parent

    @property
    def index(self):
        return int(self._index) + 1

    @property
    def core(self):
        if self._core:
            return 'M'
        else:
            return 'S'

    # TODO Consider removing any reference to OPLS, if possible
    # Otherwise, use SMIRks to find the best match
    @property
    def OPLS_type(self):
        return 'OFFT'  # stands for OpenForceField type

    # TODO Review the actual purpose of this attribute in PELE
    @property
    @none_to_zero
    def unknown(self):
        return super().unknown
