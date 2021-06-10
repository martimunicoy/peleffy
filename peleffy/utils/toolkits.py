"""
This module contains external toolkit wrappers that are required by the
main peleffy modules.
"""

import importlib
from distutils.spawn import find_executable
import tempfile
import os
import subprocess
from pathlib import Path
from copy import deepcopy

import numpy as np
from simtk import unit

from peleffy.utils import temporary_cd


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
            The path to the molecule's PDB file

        Returns
        -------
        molecule : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object
        """
        from rdkit import Chem

        return Chem.rdmolfiles.MolFromPDBFile(path, removeHs=False)

    def from_pdb_block(self, pdb_block):
        """
        It initializes an RDKit's Molecule object from a PDB block.

        Parameters
        ----------
        pdb_block : str
            The PDB block to built the molecule with

        Returns
        -------
        molecule : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object
        """
        from rdkit import Chem

        return Chem.rdmolfiles.MolFromPDBBlock(pdb_block, removeHs=False)

    def from_smiles(self, smiles):
        """
        It initializes an RDKit's Molecule object from a SMILES tag.

        Parameters
        ----------
        smiles : str
            The SMILES tag to construct the molecule structure with.

        Returns
        -------
        molecule : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object
        """
        from rdkit.Chem import AllChem as Chem

        molecule = Chem.MolFromSmiles(smiles)

        # Add hydrogens to molecule
        molecule = Chem.AddHs(molecule)

        # Generate 3D coordinates
        Chem.EmbedMolecule(molecule)

        return molecule

    def assign_connectivity_from_template(self, molecule):
        """
        It assigns the connectivity to an RDKit molecule according to the
        connectivity from an RDKit connectivity template.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        from rdkit.Chem import AllChem

        if molecule.connectivity_template is None:
            raise ValueError('A connectivity template must be previously '
                             + 'assigned to the molecule')

        rdkit_molecule = molecule.rdkit_molecule

        rdkit_molecule = AllChem.AssignBondOrdersFromTemplate(
            molecule.connectivity_template, rdkit_molecule)

        molecule._rdkit_molecule = rdkit_molecule

    def assign_stereochemistry_from_3D(self, molecule):
        """
        It assigns the stereochemistry to an RDKit molecule according to the
        3D coordinates in the PDB structure.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        from rdkit import Chem

        rdkit_molecule = molecule.rdkit_molecule
        Chem.rdmolops.AssignStereochemistryFrom3D(rdkit_molecule)

    def set_conformer(self, molecule, conformer):
        """
        It sets a new conformation to the molecule.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        conformer : an RDKit.Chem.rdchem.Conformer object
            The conformer to set to the molecule
        """

        rdkit_molecule = molecule.rdkit_molecule

        # Remove previous conformer
        rdkit_molecule.RemoveAllConformers()

        # Add current conformer
        rdkit_molecule.AddConformer(conformer, assignId=True)

    def get_residue_name(self, molecule):
        """
        It returns the name of the residue according to the RDKit molecule
        object.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        residue_name : str
            The name of the residue
        """
        rdkit_molecule = molecule.rdkit_molecule

        first_atom = list(rdkit_molecule.GetAtoms())[0]

        # Catch a None return
        try:
            residue_name = first_atom.GetPDBResidueInfo().GetResidueName()
        except AttributeError:
            residue_name = None

        return residue_name

    def get_atom_names(self, molecule):
        """
        It returns the ordered list of atom names according to the
        RDKit molecule object. In case no atom names are available
        (non-PDB source), it assignes a name to each atom considering
        the element and an index obtained from the total number of
        occurrences of each element.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        residue_name : list[str]
            The list of atom names
        """
        rdkit_molecule = molecule.rdkit_molecule

        atom_names = list()
        occurrences = dict()

        for atom in rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()

            if pdb_info is not None:
                atom_names.append(pdb_info.GetName())
            else:
                element = atom.GetSymbol()
                occurrences[element] = occurrences.get(element, 0) + 1

                atom_names.append('{:^4}'.format(str(element)
                                                 + str(occurrences[element])))

        return atom_names

    def get_dihedral(self, mol, atom1, atom2, atom3, atom4, units="radians"):
        """
        It calculates the value of the dihedral angle in the specified units
            (default radians)

        Parameters
        ----------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        atom1 : int
            Index of the first atom in the dihedral
        atom2 : int
            Index of the second atom in the dihedral
        atom3 : int
            Index of the third atom in the dihedral
        atom4 : int
            Index of the fourth atom in the dihedral
        units : str
            The units in which to calculate the angle (default is radians, can
            be radians or degrees)
        """
        from rdkit.Chem import rdMolTransforms
        if units == "degrees":
            angle = rdMolTransforms.GetDihedralDeg(mol.rdkit_molecule.GetConformer(), atom1, atom2, atom3, atom4)
        else:
            angle = rdMolTransforms.GetDihedralRad(mol.rdkit_molecule.GetConformer(), atom1, atom2, atom3, atom4)
        return angle

    def get_substruct_match(self, mol1, mol2):
        """
        It returns the atoms in mol2 that match those of mol1

        Parameters
        ----------
        mol1 : an offpele.topology.Molecule
            Molecule to match atoms from
        mol2 : an offpele.topology.Molecule
            Molecule with the atoms to match

        Returns
        -------
        mol2_atoms : tuple[int]
            The tuple of atom indices from mol2 that match those in mol1
        """
        mol2_atoms = mol1.rdkit_molecule.GetSubstructMatch(mol2.rdkit_molecule)

        return mol2_atoms

    def get_atom_degrees(self, molecule):
        """
        It returns the ordered list of atom degrees. The degree of an atom
        is defined as the number of directly-bonded neighbors. Note that
        the degree is independent of bond orders.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        atom_degrees : list[int]
            The list of atom degrees
        """
        rdkit_molecule = molecule.rdkit_molecule

        atom_degrees = list()

        for atom in rdkit_molecule.GetAtoms():
            atom_degrees.append(atom.GetDegree())

        return atom_degrees

    def get_hydrogen_parents(self, molecule):
        """
        It returns the ordered list of the element belonging to the atom
        parent for the hydrogen atoms of the molecule.

        Note that this functions sets the element to None when
        the child is not a hydrogen atom.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        atom_parents : list[int]
            The list of elements belonging to the atom parent, if the
            atom is an hydrogen. The element is set to None when the
            child is not a hydrogen atom
        """
        rdkit_molecule = molecule.rdkit_molecule

        atom_parents = list()

        for atom in rdkit_molecule.GetAtoms():
            if atom.GetSymbol() == 'H':
                bonds = atom.GetBonds()

                assert len(bonds) == 1, \
                    'Hydrogen atom should only have 1 bond'

                if bonds[0].GetBeginAtom().GetSymbol() == 'H':
                    atom_parents.append(bonds[0].GetEndAtom().GetSymbol())

                else:
                    atom_parents.append(bonds[0].GetBeginAtom().GetSymbol())

            else:
                atom_parents.append(None)

        return atom_parents

    def get_elements(self, molecule):
        """
        It returns the ordered list of elements of the molecule.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        elements : list[str]
            The list of elements belonging to supplied Molecule object
        """
        rdkit_molecule = molecule.rdkit_molecule

        elements = list()

        for atom in rdkit_molecule.GetAtoms():
            elements.append(atom.GetSymbol())

        return elements

    def to_pdb_file(self, molecule, path):
        """
        It writes the RDKit molecule to a PDB file.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        path : str
            Path to write to
        """
        from rdkit import Chem

        assert Path(path).suffix == '.pdb', 'Wrong extension'

        rdkit_molecule = molecule.rdkit_molecule

        pdb_block = Chem.rdmolfiles.MolToPDBBlock(rdkit_molecule)
        names = molecule.get_pdb_atom_names()
        tag = molecule.tag

        renamed_pdb_block = ''
        atom_counter = 0
        for line in pdb_block.split('\n'):
            if line.startswith('HETATM'):
                renamed_pdb_block += line[:12] + names[atom_counter] \
                    + ' ' + tag + line[20:] + '\n'
                atom_counter += 1
            else:
                renamed_pdb_block += line + '\n'

        with open(path, 'w') as f:
            f.write(renamed_pdb_block)

    def to_sdf_file(self, molecule, path):
        """
        It writes the RDKit molecule to an sdf file.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
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
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
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
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        rot_bonds_atom_ids : tuple[tuple[int, int]]
            The set of atom id pairs that belong to rotatable bonds
        """
        from rdkit import Chem

        rdkit_molecule = deepcopy(molecule.rdkit_molecule)

        rot_bonds_atom_ids = set([
            frozenset(atom_pair) for atom_pair in
            rdkit_molecule.GetSubstructMatches(
                Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]'))])

        # Include missing rotatable bonds for amide groups
        for atom_pair in [frozenset(atom_pair) for atom_pair in
                          rdkit_molecule.GetSubstructMatches(
                          Chem.MolFromSmarts('[$(N!@C(=O))]-&!@[!$(C(=O))&!D1&!$(*#*)]'))]:
            rot_bonds_atom_ids.add(atom_pair)

        # Remove bonds to terminal -CH3
        if molecule.exclude_terminal_rotamers:
            terminal_bonds = set([
                frozenset(atom_pair) for atom_pair in
                rdkit_molecule.GetSubstructMatches(
                    Chem.MolFromSmarts('*-&!@[$([C;H3;X4]),$([N;H2;X3]),$([N;H3;X4]),$([O;H1;X2])]'))
            ])
            rot_bonds_atom_ids = rot_bonds_atom_ids.difference(terminal_bonds)

        return list(rot_bonds_atom_ids)

    def get_coordinates(self, molecule):
        """
        It returns the 3D coordinates of all atoms in the RDKit molecule.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        coordinates : numpy.ndarray
            The array of 3D coordinates of all the atoms in the molecule
        """
        rdkit_molecule = molecule.rdkit_molecule

        conformer = rdkit_molecule.GetConformer()
        return conformer.GetPositions()

    def get_2D_representation(self, molecule):
        """
        It returns the 2D representation of the RDKit molecule.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        representation_2D : an RDKit.molecule object
            It is an RDKit molecule with an embeded 2D representation
        """
        from rdkit.Chem import AllChem

        rdkit_molecule = molecule.rdkit_molecule
        representation_2D = deepcopy(rdkit_molecule)

        AllChem.Compute2DCoords(representation_2D)
        return representation_2D

    def get_rmsd(self, molecule, molecule_2):
        """
        It returns the RMSD between two RDKit molecules.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        molecule_2 : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        rmsd_value : float
            RMSD between two RDKit molecules
        """
        from rdkit.Chem import rdMolAlign

        rmsd_value = rdMolAlign.AlignMol(molecule.rdkit_molecule,
                                         molecule_2.rdkit_molecule)
        return rmsd_value

    def draw_molecule(self, representation, atom_indexes=list(),
                      radii_dict=dict(), atom_color_dict=dict(),
                      bond_indexes=list(), bond_color_dict=dict()):
        """
        Given a molecular representation, it returns its image as SVG.

        Parameters
        ----------
        representation : an RDKit.molecule object
            It is an RDKit molecule with an embeded 2D representation

        Returns
        -------
        image : an IPython's SVG display object
            The image of the molecular representation to display
        """
        from rdkit.Chem.Draw import rdMolDraw2D

        draw = rdMolDraw2D.MolDraw2DSVG(500, 500)
        draw.SetLineWidth(4)
        rdMolDraw2D.PrepareAndDrawMolecule(draw, representation,
                                           highlightAtoms=atom_indexes,
                                           highlightAtomRadii=radii_dict,
                                           highlightAtomColors=atom_color_dict,
                                           highlightBonds=bond_indexes,
                                           highlightBondColors=bond_color_dict)
        draw.FinishDrawing()

        from IPython.display import SVG

        image = SVG(draw.GetDrawingText())

        return image


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
                'The required toolkit {} is not '.format(self.toolkit_name)
                + 'available.')

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
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
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

                self._rdkit_toolkit_wrapper.to_sdf_file(
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

        assert len(charges) == len(molecule.rdkit_molecule.GetAtoms()), \
            'Partial charge computation failed as the length of ' \
            + 'resulting partial charges does not match with the ' \
            + 'number of atoms in molecule'

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
                'The required toolkit {} is not '.format(self.toolkit_name)
                + 'available.')

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
            importlib.import_module('openff.toolkit')
            return True
        except ImportError:
            return False

    def from_rdkit(self, molecule):
        """
        It initializes an OpenForceField's Molecule object from an RDKit
        molecule.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        molecule : an openforcefield.topology.Molecule object
            The OpenForceField's Molecule
        """
        from openff.toolkit.topology.molecule import Molecule

        rdkit_molecule = molecule.rdkit_molecule
        return Molecule.from_rdkit(
            rdkit_molecule,
            allow_undefined_stereo=molecule.allow_undefined_stereo)

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
        from openff.toolkit.typing.engines.smirnoff import ForceField

        if isinstance(forcefield_name, str):
            forcefield = ForceField(forcefield_name)
        else:
            raise Exception('Invalid forcefield type')

        return forcefield

    def get_parameters_from_forcefield(self, forcefield, molecule):
        """
        It returns the parameters that are obtained with the supplied
        forcefield for a certain peleffy's molecule.

        Parameters
        ----------
        forcefield : str or an openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield from which the parameters will be obtained
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        openforcefield_parameters : dict
            The OpenFF parameters stored in a dict keyed by parameter type
        """
        from openff.toolkit.typing.engines.smirnoff import ForceField
        from openff.toolkit.topology import Topology

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

        return molecule_parameters_list[0]

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
        from openff.toolkit.typing.engines.smirnoff import ForceField

        if isinstance(forcefield, str):
            forcefield = ForceField(forcefield)
        elif isinstance(forcefield, ForceField):
            pass
        else:
            raise Exception('Invalid forcefield type')

        return forcefield.get_parameter_handler(parameter_handler_name)


class SchrodingerToolkitWrapper(ToolkitWrapper):
    """
    SchrodingerToolkitWrapper class.
    """

    _toolkit_name = 'Schrodinger Toolkit'

    def __init__(self):
        """
        It initializes a SchrodingerToolkitWrapper object.
        """
        super().__init__()

        if "SCHRODINGER" not in os.environ:
            import logging
            logging.warning("Schrodinger Toolkit requires the environment "
                            + "variable SCHRODINGER to be previously set, "
                            + "pointing to the Schrodinger's installation "
                            + "path. For more information, please, refer to "
                            + "https://martimunicoy.github.io/peleffy/installation.html#external-dependencies",
                            )

        if not self.is_available():
            raise ToolkitUnavailableException(
                'The required toolkit {} is not '.format(self.toolkit_name)
                + 'available.')

        self._rdkit_toolkit_wrapper = RDKitToolkitWrapper()

    @staticmethod
    def is_available():
        """
        Check whether the OpenForceField toolkit is installed

        Returns
        -------
        is_installed : bool
            True if OpenForceField is installed, False otherwise.
        """
        if not(RDKitToolkitWrapper.is_available()):
            return False

        if SchrodingerToolkitWrapper.path_to_ffld_server() is None:
            return False

        return True

    @staticmethod
    def path_to_ffld_server():
        FFLD_SERVER_PATH = find_executable("ffld_server")

        if FFLD_SERVER_PATH is not None:
            return FFLD_SERVER_PATH

        else:
            if "SCHRODINGER" in os.environ:
                schrodinger_root = os.environ.get('SCHRODINGER')
                return os.path.join(schrodinger_root,
                                    'utilities', 'ffld_server')

        return None

    def run_ffld_server(self, molecule):
        """
        It calls Schrodinger's ffld_server to parameterize a molecule
        with OPLS.

        .. todo ::

           * Review PlopRotTemp's atom type fixes. Should we apply them here?

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object

        Returns
        -------
        ffld_output : str
            The ffld_server output
        """

        ffld_server_exec = self.path_to_ffld_server()

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):

                self._rdkit_toolkit_wrapper.to_pdb_file(
                    molecule, tmpdir + '/molecule.pdb')

                subprocess.check_output([ffld_server_exec,
                                         "-ipdb", "molecule.pdb",
                                         "-version", "14",
                                         "-print_parameters",
                                         "-out_file", "parameters.txt"])

                with open('parameters.txt') as parameters_file:
                    return parameters_file.read()
