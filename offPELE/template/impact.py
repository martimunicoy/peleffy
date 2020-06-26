
# Global imports
from copy import deepcopy

from simtk import unit

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
        # We will work with a copy to prevent the modification of the original
        # object
        self._molecule = deepcopy(molecule)
        self._sort()

    def _sort(self):
        sorted_atoms = list()

        # Sort by core attribute and parent index
        for atom in sorted(self.molecule.atoms,
                           key=lambda a: (WritableAtom(a).core,
                                          WritableAtom(a).parent.index)):
            sorted_atoms.append(atom)

        # Define reindexer and reindex atoms
        reindexer = dict()
        for new_index, atom in enumerate(sorted_atoms):
            old_index = atom._index
            reindexer[old_index] = new_index
            atom.set_index(new_index)

        # Replace old atom list by the sorted one
        self.molecule._atoms = sorted_atoms

        # Reindex bonds, angles, propers and impropers
        for bond in self.molecule.bonds:
            bond.set_atom1_idx(reindexer[bond.atom1_idx])
            bond.set_atom2_idx(reindexer[bond.atom2_idx])
        for angle in self.molecule.angles:
            angle.set_atom1_idx(reindexer[angle.atom1_idx])
            angle.set_atom2_idx(reindexer[angle.atom2_idx])
            angle.set_atom3_idx(reindexer[angle.atom3_idx])
        for proper in self.molecule.propers:
            proper.set_atom1_idx(reindexer[proper.atom1_idx])
            proper.set_atom2_idx(reindexer[proper.atom2_idx])
            proper.set_atom3_idx(reindexer[proper.atom3_idx])
            proper.set_atom4_idx(reindexer[proper.atom4_idx])
        for improper in self.molecule.impropers:
            improper.set_atom1_idx(reindexer[improper.atom1_idx])
            improper.set_atom2_idx(reindexer[improper.atom2_idx])
            improper.set_atom3_idx(reindexer[improper.atom3_idx])
            improper.set_atom4_idx(reindexer[improper.atom4_idx])

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
        file.write('*\n')

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
            file.write('\n')
        # TODO Should we add the interactions matrix here?

    def _write_nbon(self, file):
        file.write('NBON\n')
        for atom in self.molecule.atoms:
            w_atom = WritableAtom(atom)
            # TODO an extra space is found in the IMPACT file generated by
            # PlopRotTemp, consider removing it
            file.write(' ')
            # Atom id
            file.write('{:5d}'.format(w_atom.index))
            file.write(' ')
            # Sigma
            file.write('{: 8.4f}'.format(w_atom.sigma))
            file.write(' ')
            # Epsilon
            file.write('{: 8.4f}'.format(w_atom.epsilon))
            file.write(' ')
            # Charge
            file.write('{: 10.6f}'.format(w_atom.charge))
            file.write(' ')
            # Rad. Non Polar SGB
            file.write('{: 8.4f}'.format(w_atom.born_radius))
            file.write(' ')
            # Rad. Non Polar Type
            file.write('{: 8.4f}'.format(w_atom.SASA_radius))
            file.write(' ')
            # SGB Non Polar gamma
            file.write('{: 13.9f}'.format(w_atom.nonpolar_gamma))
            file.write(' ')
            # SGB Non Polar type
            file.write('{: 13.9f}'.format(w_atom.nonpolar_alpha))
            file.write('\n')

    def _write_bond(self, file):
        file.write('BOND\n')
        for bond in self.molecule.bonds:
            w_bond = WritableBond(bond)
            idx1, idx2, spring, eq_dist = [attr[1] for attr in list(w_bond)]
            # TODO an extra space is found in the IMPACT file generated by
            # PlopRotTemp, consider removing it
            file.write(' ')
            # Atom 1 id
            file.write('{:5d}'.format(idx1))
            file.write(' ')
            # Atom 2 id
            file.write('{:5d}'.format(idx2))
            file.write(' ')
            # Spring constant (PELE works with half of the OFF's spring)
            file.write('{: 9.3f}'.format(spring / 2.0))
            file.write(' ')
            # Equilibrium distance
            file.write('{: 6.3f}\n'.format(eq_dist))

    def _write_thet(self, file):
        file.write('THET\n')
        for angle in self.molecule.angles:
            w_angle = WritableAngle(angle)
            idx1, idx2, idx3, spring, eq_angl = [attr[1] for attr in
                                                 list(w_angle)]
            # TODO an extra space is found in the IMPACT file generated by
            # PlopRotTemp, consider removing it
            file.write(' ')
            # Atom 1 id
            file.write('{:5d}'.format(idx1))
            file.write(' ')
            # Atom 2 id
            file.write('{:5d}'.format(idx2))
            file.write(' ')
            # Atom 3 id
            file.write('{:5d}'.format(idx3))
            file.write(' ')
            # Spring constant (PELE works with half of the OFF's spring)
            file.write('{: 11.5f}'.format(spring / 2.0))
            # Equilibrium angle
            file.write('{: 11.5f}\n'.format(eq_angl))

    def _write_phi(self, file):
        file.write('PHI\n')
        for proper in self.molecule.propers:
            w_proper = WritableProper(proper)
            idx1, idx2, idx3, idx4, constant, prefactor, term = \
                [attr[1] for attr in list(w_proper)]
            # Atom 1 id
            file.write('{:5d}'.format(idx1))
            file.write(' ')
            # Atom 2 id
            file.write('{:5d}'.format(idx2))
            file.write(' ')
            # Atom 3 id
            file.write('{:5d}'.format(idx3))
            file.write(' ')
            # Atom 4 id
            file.write('{:5d}'.format(idx4))
            file.write(' ')
            # Constant
            file.write('{: 9.5f}'.format(constant))
            file.write(' ')
            # Prefactor
            file.write('{: 4.1f}'.format(prefactor))
            file.write(' ')
            # Number of term
            file.write('{:3.1f}'.format(term))
            file.write('\n')

    def _write_iphi(self, file):
        file.write('IPHI\n')
        for improper in self.molecule.impropers:
            w_improper = WritableImproper(improper)
            idx1, idx2, idx3, idx4, constant, prefactor, term = \
                [attr[1] for attr in list(w_improper)]
            # TODO an extra space is found in the IMPACT file generated by
            # PlopRotTemp, consider removing it
            file.write(' ')
            # Atom 1 id
            file.write('{:5d}'.format(idx1))
            file.write(' ')
            # Atom 2 id
            file.write('{:5d}'.format(idx2))
            file.write(' ')
            # Atom 3 id
            file.write('{:5d}'.format(idx3))
            file.write(' ')
            # Atom 4 id
            file.write('{:5d}'.format(idx4))
            file.write(' ')
            # Constant
            file.write('{: 9.5f}'.format(constant))
            file.write(' ')
            # Prefactor
            file.write('{: 4.1f}'.format(prefactor))
            file.write(' ')
            # Number of term
            file.write('{:3.1f}'.format(term))
            file.write('\n')

    def _write_end(self, file):
        file.write('END\n')

    @property
    def molecule(self):
        return self._molecule


class WritableWrapper(object):
    @staticmethod
    def none_to_zero(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            if out is None:
                out = int(0)
            return out
        return function_wrapper

    @staticmethod
    def dummy_to_writable(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            out = WritableAtom(out)
            return out
        return function_wrapper

    @staticmethod
    def none_to_dummy(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            if out is None:
                out = offPELE.topology.molecule.DummyAtom(index=-1)
            return out
        return function_wrapper

    @staticmethod
    def in_angstrom(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            return out.value_in_unit(unit.angstrom)
        return function_wrapper

    @staticmethod
    def in_kcalmol(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            return out.value_in_unit(unit.kilocalorie / unit.mole)
        return function_wrapper

    @staticmethod
    def in_elementarycharge(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            return out.value_in_unit(unit.elementary_charge)
        return function_wrapper

    @staticmethod
    def in_kcal_rad2mol(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            return out.value_in_unit(unit.kilocalorie
                                     / (unit.radian**2 * unit.mole))
        return function_wrapper

    @staticmethod
    def in_deg(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            return out.value_in_unit(unit.degree)
        return function_wrapper

    @staticmethod
    def in_kcal_angstrom2mol(f):
        def function_wrapper(*args, **kwargs):
            out = f(*args, **kwargs)
            return out.value_in_unit(unit.kilocalorie
                                     / (unit.angstrom**2 * unit.mole))
        return function_wrapper


class WritableAtom(offPELE.topology.molecule.Atom, WritableWrapper):
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

    @property
    @WritableWrapper.dummy_to_writable
    @WritableWrapper.none_to_dummy
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
    @WritableWrapper.none_to_zero
    def unknown(self):
        return super().unknown

    @property
    @WritableWrapper.in_angstrom
    def sigma(self):
        return super().sigma

    @property
    @WritableWrapper.in_kcalmol
    def epsilon(self):
        return super().epsilon

    @property
    @WritableWrapper.in_elementarycharge
    def charge(self):
        return super().charge

    @property
    @WritableWrapper.none_to_zero
    def born_radius(self):
        return super().born_radius

    @property
    @WritableWrapper.in_angstrom
    def SASA_radius(self):
        return super().SASA_radius

    @property
    @WritableWrapper.none_to_zero
    def nonpolar_gamma(self):
        return super().nonpolar_gamma

    @property
    @WritableWrapper.none_to_zero
    def nonpolar_alpha(self):
        return super().nonpolar_alpha


class WritableBond(offPELE.topology.Bond, WritableWrapper):
    def __init__(self, bond):
        # We do not want to modify the original object
        bond = deepcopy(bond)

        assert isinstance(bond, (offPELE.topology.Bond,
                                 offPELE.topology.topology.Bond)), \
            'Wrong type: {}'.format(type(bond))

        super().__init__(index=bond.index, atom1_idx=bond.atom1_idx,
                         atom2_idx=bond.atom2_idx,
                         spring_constant=bond.spring_constant,
                         eq_dist=bond.eq_dist)

    @property
    def atom1_idx(self):
        return super().atom1_idx + 1

    @property
    def atom2_idx(self):
        return super().atom2_idx + 1

    @property
    @WritableWrapper.in_kcal_angstrom2mol
    def spring_constant(self):
        return super().spring_constant

    @property
    @WritableWrapper.in_angstrom
    def eq_dist(self):
        return super().eq_dist


class WritableAngle(offPELE.topology.Angle, WritableWrapper):
    def __init__(self, angle):
        # We do not want to modify the original object
        angle = deepcopy(angle)

        assert isinstance(angle, (offPELE.topology.Angle,
                                  offPELE.topology.topology.Angle)), \
            'Wrong type: {}'.format(type(angle))

        super().__init__(index=angle.index, atom1_idx=angle.atom1_idx,
                         atom2_idx=angle.atom2_idx, atom3_idx=angle.atom3_idx,
                         spring_constant=angle.spring_constant,
                         eq_angle=angle.eq_angle)

    @property
    def atom1_idx(self):
        return super().atom1_idx + 1

    @property
    def atom2_idx(self):
        return super().atom2_idx + 1

    @property
    def atom3_idx(self):
        return super().atom3_idx + 1

    @property
    @WritableWrapper.in_kcal_rad2mol
    def spring_constant(self):
        return super().spring_constant

    @property
    @WritableWrapper.in_deg
    def eq_angle(self):
        return super().eq_angle


class WritableProper(offPELE.topology.Proper, WritableWrapper):
    def __init__(self, proper):
        # We do not want to modify the original object
        proper = deepcopy(proper)

        assert isinstance(proper, (offPELE.topology.Proper,
                                   offPELE.topology.topology.Proper)), \
            'Wrong type: {}'.format(type(proper))

        super().__init__(index=proper.index, atom1_idx=proper.atom1_idx,
                         atom2_idx=proper.atom2_idx,
                         atom3_idx=proper.atom3_idx,
                         atom4_idx=proper.atom4_idx,
                         periodicity=proper.periodicity,
                         prefactor=proper.prefactor,
                         constant=proper.constant)

    @property
    def atom1_idx(self):
        return super().atom1_idx + 1

    @property
    def atom2_idx(self):
        return super().atom2_idx + 1

    @property
    def atom3_idx(self):
        return super().atom3_idx + 1

    @property
    def atom4_idx(self):
        return super().atom4_idx + 1

    @property
    @WritableWrapper.in_kcalmol
    def constant(self):
        return super().constant


class WritableImproper(offPELE.topology.Improper, WritableWrapper):
    def __init__(self, improper):
        # We do not want to modify the original object
        improper = deepcopy(improper)

        assert isinstance(improper, (offPELE.topology.Improper,
                                     offPELE.topology.topology.Improper)), \
            'Wrong type: {}'.format(type(improper))

        super().__init__(index=improper.index, atom1_idx=improper.atom1_idx,
                         atom2_idx=improper.atom2_idx,
                         atom3_idx=improper.atom3_idx,
                         atom4_idx=improper.atom4_idx,
                         periodicity=improper.periodicity,
                         prefactor=improper.prefactor,
                         constant=improper.constant)

    @property
    def atom1_idx(self):
        return super().atom1_idx + 1

    @property
    def atom2_idx(self):
        return super().atom2_idx + 1

    @property
    def atom3_idx(self):
        return super().atom3_idx + 1

    @property
    def atom4_idx(self):
        return super().atom4_idx + 1

    @property
    @WritableWrapper.in_kcalmol
    def constant(self):
        return super().constant