"""
This module contains the tests to check peleffy's molecular representations.
"""

import pytest

import tempfile
from peleffy.topology import Molecule
from peleffy.utils.toolkits import (RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper)
from peleffy.utils import get_data_file_path, temporary_cd
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
        It tests the initialization of an peleffy's Molecule representation
        from a PDB file without connectivity and a connectivity template.
        """
        # Initialize an empty Molecule object
        molecule = Molecule()
        assert molecule.connectivity_template is None, \
            'Unexpected connectivity template'

        # Initialize a Molecule from a PDB without connectivity and
        # without a connectivity template
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

        # Initialize a Molecule from a PDB without connectivity but with
        # a connectivity template
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

        # Initialize a Molecule from a PDB with connectivity and with
        # a connectivity template
        template_path = get_data_file_path(
            'ligands/BNZ.pdb')
        template = Molecule(template_path)
        ligand_path = get_data_file_path(
            'ligands/BNZ.pdb')
        molecule = Molecule(ligand_path,
                            connectivity_template=template.rdkit_molecule)

        expected_bond_ids = [(0, 1, True), (1, 2, True), (2, 3, True),
                             (3, 4, True), (4, 5, True), (0, 5, True),
                             (0, 6, False), (1, 7, False), (2, 8, False),
                             (3, 9, False), (4, 10, False), (5, 11, False)]

        for bond in molecule.rdkit_molecule.GetBonds():
            bond_id = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(),
                       bond.GetIsAromatic())
            assert bond_id in expected_bond_ids, 'Unexpected bond id ' \
                + '{}'.format(bond_id)

    def test_PDB_residue_name(self):
        """
        It tests the PDB residue name and checks for consistency with
        Molecule tag.
        """

        def check_residue_name(name):
            """Check if residue names are valid in the output PDB file"""
            with open('molecule.pdb') as f:
                for line in f:
                    if line.startswith('HETATM'):
                        assert line[17:20] == name, 'Unexpected residue name'

        ligand_path = get_data_file_path('ligands/BNZ.pdb')

        # Checking tag assignation from PDB
        molecule = Molecule(ligand_path)
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                assert molecule.tag == 'BNZ', 'Unexpected molecule tag'
                molecule.to_pdb_file('molecule.pdb')
                check_residue_name('BNZ')

        # Checking set_tag() function
        molecule = Molecule(ligand_path)
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                molecule.set_tag('TAG')
                assert molecule.tag == 'TAG', 'Unexpected molecule tag'
                molecule.to_pdb_file('molecule.pdb')
                check_residue_name('TAG')

        # Checking default tag assignment from SMILES
        molecule = Molecule(smiles='c1ccccc1')
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                assert molecule.tag == 'UNK', 'Unexpected molecule tag'
                molecule.to_pdb_file('molecule.pdb')
                check_residue_name('UNK')

        # Checking custom tag assignment from SMILES
        molecule = Molecule(smiles='c1ccccc1', tag='BEN')
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                assert molecule.tag == 'BEN', 'Unexpected molecule tag'
                molecule.to_pdb_file('molecule.pdb')
                check_residue_name('BEN')

        # Checking second custom tag assignment from SMILES
        molecule = Molecule(smiles='c1ccccc1')
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                molecule.set_tag('BNZ')
                assert molecule.tag == 'BNZ', 'Unexpected molecule tag'
                molecule.to_pdb_file('molecule.pdb')
                check_residue_name('BNZ')
