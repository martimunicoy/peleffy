"""
This module contains external toolkit wrappers that are required by the
main offpele modules.
"""

import importlib
from distutils.spawn import find_executable
import tempfile
import os
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
from simtk import unit

from offpele.utils import temporary_cd


class ToolkitUnavailableException(Exception):
    """The requested toolkit is unavailable."""
    pass


class ChargeCalculationError(Exception):
    """An external error when calculating charges"""
    pass


class ChargeMethodUnavailableError(Exception):
    """A toolkit does not support the requested partial_charge_method combination"""
    pass


class ToolkitWrapper(object):
    """
    Toolkit wrapper base class.
    """

    _is_available = None
    _toolkit_name = None

    @property
    def toolkit_name(self):
        """
        The name of the toolkit.

        Returns
        -------
        toolkit_name : str
            The name of this ToolkitWrapper object
        """
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
    """
    RDKitToolkitWrapper class.
    """

    _toolkit_name = 'RDKit Toolkit'

    def __init__(self):
        """
        It initializes a RDKitToolkitWrapper object.
        """
        super().__init__()

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {} is not '.format(self.toolkit_name)
                + 'available.')

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
        """
        It initializes an RDKit's Molecule object from a PDB file.

        Parameters
        ----------
        path : str
            The path to the molecule's PDB file.

        Returns
        -------
        molecule : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object
        """
        from rdkit import Chem

        return Chem.rdmolfiles.MolFromPDBFile(path, removeHs=False)

    def assign_stereochemistry_from_3D(self, molecule):
        """
        It assigns the stereochemistry to an RDKit molecule according to the
        3D coordinates in the PDB structure.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        """
        from rdkit import Chem

        rdkit_molecule = molecule.rdkit_molecule
        Chem.rdmolops.AssignStereochemistryFrom3D(rdkit_molecule)

    def get_residue_name(self, molecule):
        """
        It returns the name of the residue according to the RDKit molecule
        object.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object

        Returns
        -------
        residue_name : str
            The name of the residue
        """
        rdkit_molecule = molecule.rdkit_molecule

        first_atom = list(rdkit_molecule.GetAtoms())[0]
        return first_atom.GetPDBResidueInfo().GetResidueName()

    def to_sfd_file(self, molecule, path):
        """
        It writes the RDKit molecule to an sdf file.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        path : str
            Path to write to
        """
        from rdkit import Chem

        assert Path(path).suffix == '.sdf', 'Wrong extension'

        rdkit_molecule = molecule.rdkit_molecule
        with open(path, 'w') as f:
            writer = Chem.SDWriter(f)
            writer.write(rdkit_molecule)
            writer.close()

    def to_xyz_file(self, molecule, path):
        """
        It writes the RDKit molecule to an xyz file.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        path : str
            Path to write to
        """
        from rdkit import Chem

        assert Path(path).suffix == '.xyz', 'Wrong extension'

        rdkit_molecule = molecule.rdkit_molecule
        Chem.MolToXYZFile(rdkit_molecule, path)

    def get_atom_ids_with_rotatable_bonds(self, molecule):
        """
        It returns the atom ids with rotatable bonds according to the
        RDKit molecule.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object

        Returns
        -------
        rot_bonds_atom_ids : tuple[tuple[int, int]]
            The set of atom id pairs that belong to rotatable bonds
        """
        from rdkit import Chem

        rdkit_molecule = molecule.rdkit_molecule
        # Fins rotatable bond ids as in Lipinski module in RDKit
        # https://github.com/rdkit/rdkit/blob/1bf6ef3d65f5c7b06b56862b3fb9116a3839b229/rdkit/Chem/Lipinski.py#L47
        rot_bonds_atom_ids = rdkit_molecule.GetSubstructMatches(
            Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]'))

        return rot_bonds_atom_ids

    def get_coordinates(self, molecule):
        """
        It returns the 3D coordinates of all atoms in the RDKit molecule.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object

        Returns
        -------
        coordinates : numpy.ndarray
            The array of 3D coordinates of all the atoms in the molecule
        """
        rdkit_molecule = molecule.rdkit_molecule

        conformer = rdkit_molecule.GetConformer()
        return conformer.GetPositions()


class AmberToolkitWrapper(ToolkitWrapper):
    """
    AmberToolkitWrapper class.
    """

    _toolkit_name = 'Amber Toolkit'

    def __init__(self):
        """
        It initializes a AmberToolkitWrapper object.
        """
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

    def compute_partial_charges(self, molecule, method='am1bcc'):
        """
        It computes the partial charges using antechamber.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        method : str
            The name of the method to use. One of ['gasteiger', 'am1bcc'].
            If None, 'am1bcc' will be used

        Returns
        -------
        charges : simtk.unit.Quantity
            The array of partial charges

        Raises
        ------
        ChargeMethodUnavailableError if the requested charge method can not
            be handled by this toolkit
        ChargeCalculationError if the charge method is supported by this
            toolkit, but fails
        """

        SUPPORTED_CHARGE_METHODS = {'am1bcc': {'antechamber_keyword': 'bcc'},
                                    'gasteiger': {'antechamber_keyword': 'gas'}
                                    }

        if method not in SUPPORTED_CHARGE_METHODS:
            raise ChargeMethodUnavailableError(
                'partial_charge_method '
                + '{} is not available from '.format(method)
                + 'AmberToolsToolkitWrapper. Available charge methods are '
                + list(SUPPORTED_CHARGE_METHODS.keys()))

        off_molecule = molecule.off_molecule

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                net_charge = off_molecule.total_charge / \
                    unit.elementary_charge

                self._rdkit_toolkit_wrapper.to_sfd_file(
                    molecule, tmpdir + '/molecule.sdf')

                subprocess.check_output([
                    "antechamber", "-i", "molecule.sdf", "-fi", "sdf",
                    "-o", "charged.ac", "-fo", "ac", "-pf", "yes", "-dr", "n",
                    "-c",
                    SUPPORTED_CHARGE_METHODS[method]['antechamber_keyword'],
                    "-nc", str(net_charge)])
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
    """
    OpenForceFieldToolkitWrapper class.
    """

    _toolkit_name = 'OpenForceField Toolkit'

    def __init__(self):
        """
        It initializes a OpenForceFieldToolkitWrapper object.
        """
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
        """
        It initializes an OpenForceField's Molecule object from an RDKit
        molecule.

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object

        Returns
        -------
        molecule : an openforcefield.topology.Molecule object
            The OpenForceField's Molecule
        """
        from openforcefield.topology.molecule import Molecule

        rdkit_molecule = molecule.rdkit_molecule
        return Molecule.from_rdkit(rdkit_molecule)

    def get_forcefield(self, forcefield_name):
        """
        It returns the OpenForceField's object that matches with the name
        that is supplied.

        Parameters
        ----------
        forcefield_name : str
            The name of the requested forcefield

        Returns
        -------
        forcefield : an openforcefield.typing.engines.smirnoff.ForceField
                     object
            The OpenForceField's forcefield
        """
        from openforcefield.typing.engines.smirnoff import ForceField

        if isinstance(forcefield_name, str):
            forcefield = ForceField(forcefield_name)
        else:
            raise Exception('Invalid forcefield type')

        return forcefield

    def get_parameters_from_forcefield(self, forcefield, molecule):
        """
        It returns the parameters that are obtained with the supplied
        forcefield for a certain offpele's molecule.

        Parameters
        ----------
        forcefield : str or an openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield from which the parameters will be obtained
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object

        Returns
        -------
        openforcefield_parameters : an OpenForceFieldParameters object
            The OpenForceFieldParameters object
        """
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

    def get_parameter_handler_from_forcefield(self, parameter_handler_name,
                                              forcefield):
        """
        It returns a parameter handler from the forcefield based on its
        name.

        Parameters
        ----------
        parameter_handler_name : str
            The name of the parameter handler that is requested
        forcefield : an openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield from which the parameter handler will be obtained

        Returns
        -------
        parameter_handler : an openforcefield.typing.engines.smirnoff.parameters.ParameterHandler
                            object
            The ParameterHandler that was requested
        """
        from openforcefield.typing.engines.smirnoff import ForceField

        if isinstance(forcefield, str):
            forcefield = ForceField(forcefield)
        elif isinstance(forcefield, ForceField):
            pass
        else:
            raise Exception('Invalid forcefield type')

        return forcefield.get_parameter_handler(parameter_handler_name)

    class OpenForceFieldParameters(dict):
        """
        OpenForceFieldParameters class that inherits from dict.
        """

        def __init__(self, parameters_list):
            """
            It initializes an OpenForceFieldParameters object.

            parameters_list : dict
                A dictionary keyed by force type
            """
            for key, value in parameters_list.items():
                self[key] = value

        def sigmas_from_rmin_halves(func):
            """
            It converts rmin_half values to sigmas according to:
            http://ambermd.org/Questions/vdwequation.pdf
            """
            FACTOR = 0.8908987181403393  # The inverse of the sixth root of 2

            def function_wrapper(x):
                rmin_halves = func(x)

                sigmas = dict()
                for indexes, rmin_half in rmin_halves.items():
                    sigma = FACTOR * 2 * rmin_half
                    sigmas[indexes] = sigma

                return sigmas
            return function_wrapper

        def __str__(self):
            """
            It returns the readable representation string of this object.

            Returns
            -------
            output : str
                The readable representation string
            """
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
            """
            It builds the dictionary of parameters of a specific force type.

            Parameters
            ----------
            parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters of a specific force type grouped by tuples
                with the atom ids that the parameters belong to
            attribute_name : str
                The name of the attribute that is requested

            Returns
            -------
            value_by_index : dict[tuple, parameter_value]
                The parameter values that were requested grouped by the atom
                ids they belong to (arranged as a tuple)
            """
            if parameters:
                value_by_index = dict()
                for index, parameter in parameters.items():
                    value_by_index[index] = getattr(parameter, attribute_name)

                return value_by_index

        def _build_dynamic_dicts(self, parameters, attr_core_name):
            """
            It builds the dynamically the dictionaries of parameters of a
            specific force type.

            It works with the same idea as _build_dict(), however it can
            handle multiple dictionary definitions were consecutive
            parameters of the same type are found in the force type's
            parameters dictionary. It works, for example, with the multiple
            proper and improper definitions found in the OpenForceField API.
            More information in the <ProperTorsions> and <ImproperTorsions>
            sections at: https://open-forcefield-toolkit.readthedocs.io/en/0.7.0/smirnoff.html

            Parameters
            ----------
            parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters of a specific force type grouped by tuples
                with the atom ids that the parameters belong to
            attribute_name : str
                The name of the attribute that is requested

            Returns
            -------
            value_by_index : dict[tuple, parameter_value]
                The parameter values that were requested grouped by the atom
                ids they belong to (arranged as a tuple)            """
            if parameters:
                parameters_by_index = defaultdict(dict)
                all_attr_ids_found = set()
                for index, parameter in parameters.items():
                    counter = int(1)
                    attr_name = attr_core_name + str(counter)
                    while(attr_name in parameter.to_dict()):
                        all_attr_ids_found.add(counter)
                        attr_value = getattr(parameter, attr_name)
                        parameters_by_index[index][counter] = attr_value
                        counter += 1
                        attr_name = attr_core_name + str(counter)

                output_list = list()
                for attr_id in sorted(all_attr_ids_found):
                    value_by_index = dict()
                    for index in parameters.keys():
                        if attr_id in parameters_by_index[index]:
                            value_by_index[index] = \
                                parameters_by_index[index][attr_id]
                        else:
                            value_by_index[index] = None

                    output_list.append(value_by_index)

                return output_list

        # Van der Waals parameters
        def get_vdW_parameters(self):
            """
            It returns the parameters that belong to the van der Waals force
            type.

            Returns
            -------
            vdW_parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters grouped by the atom ids they belong to
                (arranged as tuples)
            """
            if 'vdW' in self:
                return self['vdW']

        def get_vdW_sigmas(self):
            """
            It gets the sigma values of the parameterized molecule.

            Returns
            -------
            sigmas : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of sigma values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_vdW_parameters()
            return self._build_dict(parameters, 'sigma')

        def get_vdW_epsilons(self):
            """
            It gets the epsilon values of the parameterized molecule.

            Returns
            -------
            epsilons : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of epsilon values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_vdW_parameters()
            return self._build_dict(parameters, 'epsilon')

        def get_vdW_rmin_halves(self):
            """
            It gets the rmin half values of the parameterized molecule.

            Returns
            -------
            rmin_halves : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of rmin half values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_vdW_parameters()
            return self._build_dict(parameters, 'rmin_half')

        @sigmas_from_rmin_halves
        def get_vdW_sigmas_from_rmin_halves(self):
            """
            It gets the rmin half values of the parameterized molecule.
            Then, a decorator converts them into sigmas.

            Returns
            -------
            sigmas : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of sigma values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            return self.get_vdW_rmin_halves()

        # Bond parameters
        def get_bond_parameters(self):
            """
            It returns the parameters that belong to the bonding force type.

            Returns
            -------
            bond_parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters grouped by the atom ids they belong to
                (arranged as tuples)
            """
            if 'Bonds' in self:
                return self['Bonds']

        def get_bond_lengths(self):
            """
            It gets the bond length values of the parameterized molecule.

            Returns
            -------
            bond_lengths : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of bond length values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_bond_parameters()
            return self._build_dict(parameters, 'length')

        def get_bond_ks(self):
            """
            It gets the bond k values of the parameterized molecule.

            Returns
            -------
            bond_ks : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of bond k values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_bond_parameters()
            return self._build_dict(parameters, 'k')

        # Angle parameters
        def get_angle_parameters(self):
            """
            It returns the parameters that belong to the angular force type.

            Returns
            -------
            angle_parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters grouped by the atom ids they belong to
                (arranged as tuples)
            """
            if 'Angles' in self:
                return self['Angles']

        def get_angle_angles(self):
            """
            It gets the angle values of the parameterized molecule.

            Returns
            -------
            angles : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of angle values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_angle_parameters()
            return self._build_dict(parameters, 'angle')

        def get_angle_ks(self):
            """
            It gets the angle k values of the parameterized molecule.

            Returns
            -------
            angle_ks : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of angle k values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_angle_parameters()
            return self._build_dict(parameters, 'k')

        # Dihedral parameters
        def get_dihedral_parameters(self):
            """
            It returns the parameters that belong to the proper dihedrals
            force type.

            Returns
            -------
            proper_parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters grouped by the atom ids they belong to
                (arranged as tuples)
            """
            if 'ProperTorsions' in self:
                return self['ProperTorsions']

        def get_dihedral_periodicities(self):
            """
            It gets the dihedral periodicity values of the parameterized
            molecule.

            Returns
            -------
            dihedral_periodicities : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of dihedral periodicity values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_dihedral_parameters()
            return self._build_dynamic_dicts(parameters, 'periodicity')

        def get_dihedral_phases(self):
            """
            It gets the dihedral phase values of the parameterized
            molecule.

            Returns
            -------
            dihedral_phases : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of dihedral phase values grouped
                by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_dihedral_parameters()
            return self._build_dynamic_dicts(parameters, 'phase')

        def get_dihedral_ks(self):
            """
            It gets the dihedral periodicity values of the parameterized
            molecule.

            Returns
            -------
            dihedral_periodicities : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of dihedral periodicity values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_dihedral_parameters()
            return self._build_dynamic_dicts(parameters, 'k')

        def get_dihedral_idivfs(self):
            """
            It gets the dihedral idivf values of the parameterized
            molecule.

            Returns
            -------
            dihedral_idivfs : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of dihedral idivf values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_dihedral_parameters()
            return self._build_dynamic_dicts(parameters, 'idivf')

        # Improper parameters
        def get_improper_parameters(self):
            """
            It returns the parameters that belong to the improper dihedrals
            force type.

            Returns
            -------
           improper_parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters grouped by the atom ids they belong to
                (arranged as tuples)
            """
            if 'ImproperTorsions' in self:
                return self['ImproperTorsions']

        def get_improper_periodicities(self):
            """
            It gets the improper periodicity values of the parameterized
            molecule.

            Returns
            -------
            improper_periodicities : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of improper periodicity values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_improper_parameters()
            return self._build_dynamic_dicts(parameters, 'periodicity')

        def get_improper_phases(self):
            """
            It gets the improper phase values of the parameterized
            molecule.

            Returns
            -------
            improper_phases : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of improper phase values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_improper_parameters()
            return self._build_dynamic_dicts(parameters, 'phase')

        def get_improper_ks(self):
            """
            It gets the improper k values of the parameterized
            molecule.

            Returns
            -------
            improper_ks : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of improper k values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_improper_parameters()
            return self._build_dynamic_dicts(parameters, 'k')

        def get_improper_idivfs(self):
            """
            It gets the improper idivf values of the parameterized
            molecule.

            Returns
            -------
            improper_idivfs : list[dict[tuple[int], simtk.unit.Quantity]]
                The list of dictionaries of improper idivf values
                grouped by the atom ids they belong to (arranged as tuples)
            """
            parameters = self.get_improper_parameters()
            return self._build_dynamic_dicts(parameters, 'idivf')

        # GBSA solvent parameters
        def get_GBSA_parameters(self):
            """
            It returns the parameters that belong to the GBSA force type.

            Returns
            -------
           GSBA_parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters grouped by the atom ids they belong to
                (arranged as tuples)
            """
            if 'GBSA' in self:
                return self['GBSA']

        def get_GBSA_radii(self):
            """
            It gets the GBSA radius values of the parameterized molecule.

            Returns
            -------
            GBSA_radii : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of GBSA radius values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_GBSA_parameters()
            return self._build_dict(parameters, 'radius')

        def get_GBSA_scales(self):
            """
            It gets the GBSA scale values of the parameterized molecule.

            Returns
            -------
            GBSA_scales : dict[tuple[int], simtk.unit.Quantity]
                The dictionary of GBSA scale values grouped by the atom ids
                they belong to (arranged as tuples)
            """
            parameters = self.get_GBSA_parameters()
            return self._build_dict(parameters, 'scale')
