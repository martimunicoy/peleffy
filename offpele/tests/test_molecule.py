"""
This module contains the tests to check offpele's molecular representations.
"""

import pytest

from offpele.topology import Molecule
from offpele.utils.toolkits import (RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper)
from offpele.utils import get_data_file_path
from .utils import SET_OF_LIGAND_PATHS


class TestMolecule(object):
    """
    It wraps all tests that involve the Molecule class.
    """

    def test_bad_init_parameterization(self):
        """
        It checks that a call to the parameterize function with a Molecule
        unsuccessfully initialized raises an Exception.
        """
        FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'
        LIGAND_PATH = SET_OF_LIGAND_PATHS[0]
        ligand_path = get_data_file_path(LIGAND_PATH)

        molecule = Molecule()
        with pytest.raises(Exception):
            molecule.parameterize(FORCEFIELD_NAME)

        rdkit_toolkit = RDKitToolkitWrapper()
        molecule._rdkit_molecule = rdkit_toolkit.from_pdb(ligand_path)
        molecule._off_molecule = None

        with pytest.raises(Exception):
            molecule.parameterize(FORCEFIELD_NAME)

        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        molecule._off_molecule = openforcefield_toolkit.from_rdkit(molecule)
        molecule._rdkit_molecule = None

        with pytest.raises(Exception):
            molecule.parameterize(FORCEFIELD_NAME)

    def test_good_init_parameterization(self):
        """
        It checks that a call to the parameterize function with a Molecule
        successfully initialized does not raise any Exception.
        """
        FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'
        LIGAND_PATH = SET_OF_LIGAND_PATHS[0]
        ligand_path = get_data_file_path(LIGAND_PATH)

        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME)

    def test_molecule_name_assignment(self):
        """
        It tests the molecule name assignment.
        """
        # Look for an empty name when dummy Molecule is loaded
        molecule = Molecule()
        assert molecule.name == '', 'Unexpected atom name'

        # Look for the PDB name when a Molecule is loaded from a PDB file
        ligand_path = get_data_file_path('ligands/BNZ.pdb')
        molecule = Molecule(ligand_path)
        assert molecule.name == 'BNZ', 'Unexpected atom name'

        # Look for benzene name when a Molecule is loaded from a PDB file
        # with a custom name
        ligand_path = get_data_file_path('ligands/BNZ.pdb')
        molecule = Molecule(ligand_path, name='benzene')
        assert molecule.name == 'benzene', 'Unexpected atom name'

        # Look for the SMILES name when a Molecule is loaded from a SMILES tag
        molecule = Molecule(smiles='c1ccccc1')
        assert molecule.name == 'c1ccccc1', 'Unexpected atom name'

        # Look for benzene name when a Molecule is loaded from a SMILES tag
        # with a custom name
        molecule = Molecule(smiles='c1ccccc1', name='benzene')
        assert molecule.name == 'benzene', 'Unexpected atom name'

    def test_molecule_tag_assignment(self):
        """
        It tests the molecule tag assignment.
        """
        # Look for UNK tag when dummy Molecule is loaded
        molecule = Molecule()
        assert molecule.tag == 'UNK', 'Unexpected atom tag'

        # Look for the PDB residue name as a tag when a Molecule is loaded
        # from a PDB file
        ligand_path = get_data_file_path('ligands/BNZ.pdb')
        molecule = Molecule(ligand_path)
        assert molecule.tag == 'BNZ', 'Unexpected atom tag'

        # Look for BEN tag when a Molecule is loaded from a PDB file with
        # a custom name
        ligand_path = get_data_file_path('ligands/BNZ.pdb')
        molecule = Molecule(ligand_path, tag='BEN')
        assert molecule.tag == 'BEN', 'Unexpected atom tag'

        # Look for UNK tag when a Molecule is loaded from a SMILES tag
        molecule = Molecule(smiles='c1ccccc1')
        assert molecule.tag == 'UNK', 'Unexpected atom tag'

        # Look for BNZ tag when a Molecule is loaded from a SMILES tag with
        # a custom tag
        molecule = Molecule(smiles='c1ccccc1', tag='BNZ')
        assert molecule.tag == 'BNZ', 'Unexpected atom tag'

    def test_PDB_connectivity_template(self):
        """
        It tests the initialization of an offpele's Molecule representation
        from a PDB file without connectivity and a connectivity template.
        """
        # Initialize an empty Molecule object
        molecule = Molecule()
        assert molecule.connectivity_template is None, \
            'Unexpected connectivity template'

        # Initialize a Molecule from a PDB without connectivity
        ligand_path = get_data_file_path(
            'ligands/BNZ_without_connectivity.pdb')
        molecule = Molecule(ligand_path)

        expected_bond_ids = [(1, 0, False), (2, 1, False), (3, 2, False),
                             (4, 3, False), (5, 4, False), (5, 0, False),
                             (6, 0, False), (7, 1, False), (8, 2, False),
                             (9, 3, False), (10, 4, False), (11, 5, False)]

        for bond in molecule.rdkit_molecule.GetBonds():
            bond_id = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                       bond.GetIsAromatic())
            assert bond_id in expected_bond_ids, 'Unexpected bond id ' \
                + '{}'.format(bond_id)

        # Initialize a Molecule from a PDB without connectivity
        template_path = get_data_file_path(
            'ligands/BNZ.pdb')
        template = Molecule(template_path)
        ligand_path = get_data_file_path(
            'ligands/BNZ_without_connectivity.pdb')
        molecule = Molecule(ligand_path,
                            connectivity_template=template.rdkit_molecule)

        expected_bond_ids = [(1, 0, True), (2, 1, True), (3, 2, True),
                             (4, 3, True), (5, 4, True), (5, 0, True),
                             (6, 0, False), (7, 1, False), (8, 2, False),
                             (9, 3, False), (10, 4, False), (11, 5, False)]

        for bond in molecule.rdkit_molecule.GetBonds():
            bond_id = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                       bond.GetIsAromatic())
            assert bond_id in expected_bond_ids, 'Unexpected bond id ' \
                + '{}'.format(bond_id)
