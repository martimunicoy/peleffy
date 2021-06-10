"""
This module contains the tests to check peleffy's molecular representations.
"""

import pytest

import tempfile

from peleffy.topology import Molecule
from peleffy.utils import get_data_file_path, temporary_cd


class TestMolecule(object):
    """
    It wraps all tests that involve the Molecule class.
    """

    def test_pdb_initialization(self):
        """
        It checks the initialization from a PDB file.
        """
        ligand_path = get_data_file_path('ligands/ethylene.pdb')

        molecule = Molecule(ligand_path)

        # Save it
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                molecule.to_pdb_file('molecule.pdb')

    def test_smiles_initialization(self):
        """
        It checks the initialization from a SMILES tag.
        """
        molecule = Molecule(smiles='c1ccccc1')

        # Save it
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                molecule.to_pdb_file('molecule.pdb')

    def test_molecule_name_assignment(self):
        """
        It tests the molecule name assignment.
        """
        # Look for an empty name when dummy Molecule is loaded
        molecule = Molecule()
        assert molecule.name == '', 'Unexpected atom name'

        # Look for the PDB name when a Molecule is loaded from a PDB file
        ligand_path = get_data_file_path('ligands/benzene.pdb')
        molecule = Molecule(ligand_path)
        assert molecule.name == 'benzene', 'Unexpected atom name'

        # Look for benzene name when a Molecule is loaded from a PDB file
        # with a custom name
        ligand_path = get_data_file_path('ligands/benzene.pdb')
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
        ligand_path = get_data_file_path('ligands/benzene.pdb')
        molecule = Molecule(ligand_path)
        assert molecule.tag == 'BNZ', 'Unexpected atom tag'

        # Look for BEN tag when a Molecule is loaded from a PDB file with
        # a custom name
        ligand_path = get_data_file_path('ligands/benzene.pdb')
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
            'ligands/benzene_without_connectivity.pdb')
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
            'ligands/benzene.pdb')
        template = Molecule(template_path)
        ligand_path = get_data_file_path(
            'ligands/benzene_without_connectivity.pdb')
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
            'ligands/benzene.pdb')
        template = Molecule(template_path)
        ligand_path = get_data_file_path(
            'ligands/benzene.pdb')
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

        ligand_path = get_data_file_path('ligands/benzene.pdb')

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

    def test_pdb_checkup(self):
        """It tests the safety check function for PDB files."""

        LIGAND_GOOD = get_data_file_path('ligands/ethylene.pdb')
        LIGAND_ERROR1 = get_data_file_path('tests/ethylene_error1.pdb')
        LIGAND_ERROR2 = get_data_file_path('tests/ethylene_error2.pdb')
        LIGAND_ERROR3 = get_data_file_path('tests/ethylene_error3.pdb')
        LIGAND_ERROR4 = get_data_file_path('tests/ethylene_error4.pdb')

        # This should work without any complain
        _ = Molecule(LIGAND_GOOD)

        # All atom names need to be unique
        with pytest.raises(Exception):
            _ = Molecule(LIGAND_ERROR1)

        # All residue ids must match
        with pytest.raises(Exception):
            _ = Molecule(LIGAND_ERROR2)

        # All residue names must match
        with pytest.raises(Exception):
            _ = Molecule(LIGAND_ERROR3)

        # Check warning message in the logger when connectivity is missing
        import io
        from peleffy.utils import Logger
        import logging
        from importlib import reload
        logging.shutdown()
        reload(logging)

        log = Logger()
        log.set_level('WARNING')

        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            _ = Molecule(LIGAND_ERROR4)

            output = buf.getvalue()

            assert output == "Warning: input PDB has no information " \
                + "about the connectivity and this could result in " \
                + "an unexpected bond assignment\n"

    def test_undefined_stereo(self):
        """
        It checks the behaviour when ignoring the stereochemistry
        in the Molecule initialization.
        """
        from openff.toolkit.utils.toolkits import UndefinedStereochemistryError
        from peleffy.forcefield import OpenForceField

        # This should crash due to an undefined stereochemistry error
        with pytest.raises(UndefinedStereochemistryError):
            mol = Molecule(smiles='CN(C)CCC=C1c2ccccc2CCc3c1cccc3')

        # This now should work
        mol = Molecule(smiles='CN(C)CCC=C1c2ccccc2CCc3c1cccc3',
                       allow_undefined_stereo=True)

        # And we can parameterize it
        ff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        ff.parameterize(mol, charge_method='gasteiger')

        # Save it
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                mol.to_pdb_file('molecule.pdb')

    def test_from_rdkit(self):
        """
        It checks the initialization of a peleffy Molecule from an RDKit
        molecular representation.
        """
        from rdkit import Chem

        pdb_path = get_data_file_path('ligands/malonate.pdb')
        rdkit_molecule = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        molecule = Molecule.from_rdkit(rdkit_molecule)

        assert molecule._rdkit_molecule is not None, \
            'Unexpected molecule representation found, it is not initialized'
        assert molecule._off_molecule is not None, \
            'Unexpected molecule representation found, it is not initialized'

        ref_pdb_atom_names = [' O1 ', ' C1 ', ' O2 ', ' C2 ', ' C3 ',
                              ' O3 ', ' O4 ', ' H1 ', ' H2 ', ' H3 ']
        pdb_atom_names = molecule.get_pdb_atom_names()

        assert pdb_atom_names == ref_pdb_atom_names, \
            'Unexpected PDB atom names found in the resulting Molecule ' \
            + 'representation'

        assert molecule.graph is not None, \
            'Molecule\' graph should be initialized'

    def test_from_openff(self):
        """
        It checks the initialization of a peleffy Molecule from an OpenFF
        molecular representation.
        """
        from openff.toolkit.topology import Molecule as OpenFFMolecule

        openff_molecule = OpenFFMolecule.from_smiles('C(C(=O)[O-])C(=O)[OH]')

        molecule = Molecule.from_openff(openff_molecule)

        assert molecule._rdkit_molecule is not None, \
            'Unexpected molecule representation found, it is not initialized'
        assert molecule._off_molecule is not None, \
            'Unexpected molecule representation found, it is not initialized'

        ref_pdb_atom_names = [' C1 ', ' C2 ', ' O1 ', ' O2 ', ' C3 ',
                              ' O3 ', ' O4 ', ' H1 ', ' H2 ', ' H3 ']
        pdb_atom_names = molecule.get_pdb_atom_names()

        assert pdb_atom_names == ref_pdb_atom_names, \
            'Unexpected PDB atom names found in the resulting Molecule ' \
            + 'representation'

        assert molecule.graph is not None, \
            'Molecule\' graph should be initialized'

    def test_pdb_fixer(self):
        """
        It checks the PDB fixer prior parsing a PDB input file for
        a peleffy Molecule.
        """
        from .utils import compare_blocks

        # Check default
        molecule = Molecule()

        assert molecule.fix_pdb is True, \
            'Unexpected default settings for the PDB fixer'

        # Activate fixer
        molecule = Molecule(fix_pdb=True)

        ref_path = get_data_file_path('tests/ligSUV_fixed.pdb')
        path1 = get_data_file_path('tests/ligSUV_no_elements1.pdb')
        path2 = get_data_file_path('tests/ligSUV_no_elements2.pdb')

        ref_pdb_block = molecule._read_and_fix_pdb(ref_path)
        pdb_block1 = molecule._read_and_fix_pdb(path1)
        pdb_block2 = molecule._read_and_fix_pdb(path2)

        compare_blocks(ref_pdb_block, pdb_block1, (76, 78))
        compare_blocks(ref_pdb_block, pdb_block2, (76, 78))

        # Deactivate fixer
        molecule = Molecule(fix_pdb=False)

        ref_path = get_data_file_path('tests/ligSUV_fixed.pdb')
        path1 = get_data_file_path('tests/ligSUV_no_elements1.pdb')
        path2 = get_data_file_path('tests/ligSUV_no_elements2.pdb')

        ref_pdb_block = molecule._read_and_fix_pdb(ref_path)
        pdb_block1 = molecule._read_and_fix_pdb(path1)
        pdb_block2 = molecule._read_and_fix_pdb(path2)

        with pytest.raises(AssertionError):
            compare_blocks(ref_pdb_block, pdb_block1, (76, 78))

        with pytest.raises(AssertionError):
            compare_blocks(ref_pdb_block, pdb_block2, (76, 78))

    def test_pdb_fixer_logger_messages(self):
        """It checks the logger messages of the PDB fixer."""

        from peleffy.utils import Logger
        import io

        molecule = Molecule(fix_pdb=True)

        # Check logger messages
        path3 = get_data_file_path('tests/ligSUV_no_elements3.pdb')

        import logging

        # Force a hard reset of logging library and the logger it manages
        from importlib import reload
        logging.shutdown()
        reload(logging)

        # Initiate logger
        log = Logger()

        # Try the default level (INFO)
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            _ = molecule._read_and_fix_pdb(path3)

            # Get string from buffer
            output = buf.getvalue()

            assert output == "Warning: input PDB has no information " \
                + "about atom elements and they were inferred from " \
                + "atom names. " \
                + "Please, verify that the resulting elements are " \
                + "correct\n" \
                + "Error: PDB could not be fixed\n"

    def test_molecule_display(self):
        """
        It checks the visual representation of the molecule in a
        Jupyter notebook.
        """

        from IPython.display import display

        molecule = Molecule(smiles='c1ccccc1')

        # This should not raise any Exception
        display(molecule)
