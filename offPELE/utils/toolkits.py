import importlib
from distutils.spawn import find_executable
import tempfile
import os
import subprocess

import numpy as np
from simtk import unit

from offPELE.utils import temporary_cd


# Exceptions
class ToolkitUnavailableException(Exception):
    """The requested toolkit is unavailable."""
    pass


class ChargeCalculationError(Exception):
    """An external error when calculating charges"""
    pass


class ToolkitWrapper(object):
    """
    Toolkit wrapper base class.
    """
    _is_available = None
    _toolkit_name = None

    @property
    def toolkit_name(self):
        return self._toolkit_name

    @staticmethod
    def is_available():
        """
        Check whether the corresponding toolkit can be imported
        Returns
        -------
        is_installed : bool
            True if corresponding toolkit is installed, False otherwise.
        """
        return NotImplementedError


class RDKitToolkitWrapper(ToolkitWrapper):
    _toolkit_name = 'RDKit Toolkit'

    def __init__(self):
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {self.toolkit_name} is not '
                'available.')

    @staticmethod
    def is_available():
        """
        Check whether the RDKit toolkit can be imported
        Returns
        -------
        is_installed : bool
            True if RDKit is installed, False otherwise.
        """
        try:
            importlib.import_module('rdkit', 'Chem')
            return True
        except ImportError:
            return False

    def from_pdb(self, path):
        from rdkit import Chem

        return Chem.rdmolfiles.MolFromPDBFile(path, removeHs=False)

    def assign_stereochemistry_from_3D(self, molecule):
        from rdkit import Chem

        rdkit_molecule = molecule.rdkit_molecule
        Chem.rdmolops.AssignStereochemistryFrom3D(rdkit_molecule)

    def to_sfd_file(self, molecule, path):
        from rdkit import Chem

        rdkit_molecule = molecule.rdkit_molecule
        with open('molecule.sdf', 'w') as f:
            writer = Chem.SDWriter(f)
            writer.write(rdkit_molecule)
            writer.close()


class AmberToolkitWrapper(ToolkitWrapper):
    _toolkit_name = 'Amber Toolkit'

    def __init__(self):
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {self.toolkit_name} is not '
                'available.')

        self._rdkit_toolkit_wrapper = RDKitToolkitWrapper()

    @staticmethod
    def is_available():
        """
        Check whether the AmberTools toolkit is installed
        Returns
        -------
        is_installed : bool
            True if AmberTools is installed, False otherwise.
        """
        ANTECHAMBER_PATH = find_executable("antechamber")
        if ANTECHAMBER_PATH is None:
            return False
        if not(RDKitToolkitWrapper.is_available()):
            return False
        return True

    def compute_partial_charges_am1bcc(self, molecule):
        off_molecule = molecule.off_molecule

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = off_molecule.total_charge / \
                    unit.elementary_charge

                self._rdkit_toolkit_wrapper.to_sfd_file(molecule, tmpdir)

                subprocess.check_output([
                    "antechamber", "-i", "molecule.sdf", "-fi", "sdf",
                    "-o", "charged.ac", "-fo", "ac", "-pf", "yes", "-dr", "n",
                    "-c", "bcc", "-nc", str(net_charge)])
                # Write out just charges
                subprocess.check_output([
                    "antechamber", "-dr", "n", "-i", "charged.ac", "-fi", "ac",
                    "-o", "charged2.ac", "-fo", "ac", "-c", "wc",
                    "-cf", "charges.txt", "-pf", "yes"])

                if not os.path.exists('charges.txt'):
                    # TODO: copy files into local directory to aid debugging?
                    raise ChargeCalculationError(
                        "Antechamber/sqm partial charge calculation failed on "
                        "molecule {} (SMILES {})".format(
                            off_molecule.name, off_molecule.to_smiles()))

                # Read the charges
                with open('charges.txt', 'r') as infile:
                    contents = infile.read()

                text_charges = contents.split()
                charges = np.zeros([off_molecule.n_atoms], np.float64)
                for index, token in enumerate(text_charges):
                    charges[index] = float(token)

        charges = unit.Quantity(charges, unit.elementary_charge)

        return charges


class OpenForceFieldToolkitWrapper(ToolkitWrapper):
    _toolkit_name = 'OpenForceField Toolkit'

    def __init__(self):
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {self.toolkit_name} is not '
                'available.')

    @staticmethod
    def is_available():
        """
        Check whether the OpenForceField toolkit is installed
        Returns
        -------
        is_installed : bool
            True if OpenForceField is installed, False otherwise.
        """
        try:
            importlib.import_module('openforcefield')
            return True
        except ImportError:
            return False

    def from_rdkit(self, molecule):
        from openforcefield.topology.molecule import Molecule

        rdkit_molecule = molecule.rdkit_molecule
        return Molecule.from_rdkit(rdkit_molecule)

    def get_parameters_from_forcefield(self, forcefield, molecule):
        from openforcefield.typing.engines.smirnoff import ForceField
        from openforcefield.topology import Topology

        off_molecule = molecule.off_molecule
        topology = Topology.from_molecules([off_molecule])

        if isinstance(forcefield, str):
            forcefield = ForceField(forcefield)
        elif isinstance(forcefield, ForceField):
            pass
        else:
            raise Exception('Invalid forcefield type')

        molecule_parameters_list = forcefield.label_molecules(topology)

        assert len(molecule_parameters_list) == 1, 'A single molecule is ' \
            'expected'
        return self.OpenForceFieldParameters(molecule_parameters_list[0])

    class OpenForceFieldParameters(dict):
        def __init__(self, parameters_list):
            for key, value in parameters_list.items():
                self[key] = value

        def __str__(self):
            output = ''
            for force_tag, force_dict in self.items():
                output += f"\n{force_tag}:\n"
                for (atom_indices, parameter) in force_dict.items():
                    atomstr = ''
                    for idx in atom_indices:
                        atomstr += '%3s' % idx
                    output += " - atoms: %s  parameter_id: %s  smirks %s\n" % \
                        (atomstr, parameter.id, parameter.smirks)
            return output

        def _build_dict(self, parameters, attribute_name):
            output_dict = dict()
            for index, parameter in parameters.items():
                try:
                    output_dict[index] = getattr(parameter, attribute_name)
                except IndexError:
                    output_dict[index] = None

            return output_dict

        # Van der Waals parameters
        def get_vdW_parameters(self):
            if 'vdW' in self:
                return self['vdW']

        def get_vdW_sigmas(self):
            parameters = self.get_vdW_parameters()
            return self._build_dict(parameters, 'sigma')

        def get_vdW_epsilons(self):
            parameters = self.get_vdW_parameters()
            return self._build_dict(parameters, 'epsilon')

        def get_vdW_rmin_halves(self):
            parameters = self.get_vdW_parameters()
            return self._build_dict(parameters, 'rmin_half')

        # Bond parameters
        def get_bond_parameters(self):
            if 'Bonds' in self:
                return self['Bonds']

        def get_bond_lengths(self):
            parameters = self.get_bond_parameters()
            return self._build_dict(parameters, 'length')

        def get_bond_ks(self):
            parameters = self.get_bond_parameters()
            return self._build_dict(parameters, 'k')

        # Angle parameters
        def get_angle_parameters(self):
            if 'Angles' in self:
                return self['Angles']

        def get_angle_angles(self):
            parameters = self.get_angle_parameters()
            return self._build_dict(parameters, 'angle')

        def get_angle_ks(self):
            parameters = self.get_angle_parameters()
            return self._build_dict(parameters, 'k')

        # Dihedral parameters
        def get_dihedral_parameters(self):
            if 'ProperTorsions' in self:
                return self['ProperTorsions']

        def get_dihedral_periodicity1s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'periodicity1')

        def get_dihedral_periodicity2s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'periodicity2')

        def get_dihedral_periodicity3s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'periodicity3')

        def get_dihedral_phase1s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'phase1')

        def get_dihedral_phase2s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'phase2')

        def get_dihedral_phase3s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'phase3')

        def get_dihedral_k1s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'k1')

        def get_dihedral_k2s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'k2')

        def get_dihedral_k3s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'k3')

        def get_dihedral_idivf1s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'idivf1')

        def get_dihedral_idivf2s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'idivf2')

        def get_dihedral_idivf3s(self):
            parameters = self.get_dihedral_parameters()
            return self._build_dict(parameters, 'idivf3')

        # Improper parameters
        def get_improper_parameters(self):
            if 'ImproperTorsions' in self:
                return self['ImproperTorsions']

        def get_improper_periodicity1s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'periodicity1')

        def get_improper_periodicity2s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'periodicity2')

        def get_improper_periodicity3s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'periodicity3')

        def get_improper_phase1s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'phase1')

        def get_improper_phase2s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'phase2')

        def get_improper_phase3s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'phase3')

        def get_improper_k1s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'k1')

        def get_improper_k2s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'k2')

        def get_improper_k3s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'k3')

        def get_improper_idiv1s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'idiv1')

        def get_improper_idiv2s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'idiv2')

        def get_improper_idiv3s(self):
            parameters = self.get_improper_parameters()
            return self._build_dict(parameters, 'idiv3')
