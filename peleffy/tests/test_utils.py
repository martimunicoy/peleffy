"""
This module contains the tests to check some handy classes and functions
of peleffy.
"""


import pytest

import io
import os
import tempfile

from peleffy.utils import Logger
from peleffy.topology import Molecule


class TestLogger(object):
    def test_logger_levels(self):
        """
        It checks the correct behaviour of the different log levels.
        """
        def push_messages(log):
            """Pull some messages at different levels."""
            log.debug('Debug message')
            log.info('Info message')
            log.warning('Warn message')
            log.error('Error message')
            log.critical('Critical message')

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

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Info message\nWarn message\n' \
                + 'Error message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try DEBUG level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try DEBUG level
            log.set_level('DEBUG')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Debug message\nInfo message\n'\
                + 'Warn message\nError message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try INFO level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try INFO level
            log.set_level('INFO')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Info message\nWarn message\n' \
                + 'Error message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try WARNING level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try WARNING level
            log.set_level('WARNING')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Warn message\nError message\n' \
                + 'Critical message\n', \
                'Unexpected logger message at standard output'

        # Try ERROR level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try ERROR level
            log.set_level('ERROR')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Error message\nCritical message\n', \
                'Unexpected logger message at standard output'

        # Try CRITICAL level
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            # Try CRITICAL level
            log.set_level('CRITICAL')

            # Push messages
            push_messages(log)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Critical message\n', \
                'Unexpected logger message at standard output'


class TestOutputPathHandler(object):
    """
    It contains all the tests to validate the OuputPathHandler class.
    """

    def test_non_datalocal_paths(self):
        """
        It tests the non-datalocal paths assignment.
        """
        from peleffy.utils import OutputPathHandler
        from peleffy.forcefield import OpenForceField

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')

        # Load force field
        openff = OpenForceField('openff_unconstrained-1.2.1.offxml')

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, openff,
                                           as_datalocal=False)

        # Validate output paths
        assert output_handler.get_rotamer_library_path() == \
            './{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path() == \
            './{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path() == \
            './ligandParams.txt', \
            'Unexpected default solvent parameters path'
        assert output_handler.get_conformation_library_path() == \
            './{}.conformation'.format(tag.upper()), \
            'Unexpected default conformation library path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, openff, as_datalocal=False,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path() == \
                os.path.join(tmpdir, 'output',
                             '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path() == \
                os.path.join(tmpdir, 'output', '{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path() == \
                os.path.join(tmpdir, 'output', 'ligandParams.txt'), \
                'Unexpected default solvent parameters path'
            assert output_handler.get_conformation_library_path() == \
                os.path.join(tmpdir, 'output',
                             '{}.conformation'.format(tag.upper())), \
                'Unexpected default conformation library path'

    def test_datalocal_paths_for_openff(self):
        """It tests the datalocal paths assignment for OpenFF."""
        from peleffy.utils import OutputPathHandler
        from peleffy.forcefield import OpenForceField

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')

        # Load force field
        openff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, openff,
                                           as_datalocal=True)

        # Validate output paths
        assert output_handler.get_rotamer_library_path(
            create_missing_folders=False) == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path(
            create_missing_folders=False) == \
            './DataLocal/Templates/OpenFF/Parsley/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path(
            create_missing_folders=False) == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'
        assert output_handler.get_conformation_library_path(
            create_missing_folders=False) == \
            './DataLocal/Conformations/' \
            + '{}.conformation'.format(tag.upper()), \
            'Unexpected default conformation library path'

        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, openff, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            # Validate output paths
            assert output_handler.get_rotamer_library_path(
                create_missing_folders=False) == \
                tmpdir + '/output/DataLocal/LigandRotamerLibs/' \
                + '{}.rot.assign'.format(tag.upper()), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path(
                create_missing_folders=False) == \
                tmpdir + '/output/DataLocal/Templates/OpenFF/Parsley/' \
                + '{}z'.format(tag.lower()), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path(
                create_missing_folders=False) == \
                tmpdir + '/output/DataLocal/OBC/ligandParams.txt', \
                'Unexpected default solvent parameters path'
            assert output_handler.get_conformation_library_path(
                create_missing_folders=False) == \
                tmpdir + '/output/DataLocal/Conformations/' \
                + '{}.conformation'.format(tag.upper()), \
                'Unexpected default conformation library path'

    def test_datalocal_paths_for_opls(self):
        """It tests the datalocal paths assignment for OPLS2005."""
        from peleffy.utils import OutputPathHandler
        from peleffy.forcefield import OPLS2005ForceField

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')

        # Load force field
        opls2005 = OPLS2005ForceField()

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, opls2005,
                                           as_datalocal=True)

        # Validate output paths
        assert output_handler.get_rotamer_library_path(
            create_missing_folders=False) == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path(
            create_missing_folders=False) == \
            './DataLocal/Templates/OPLS2005/HeteroAtoms/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path(
            create_missing_folders=False) == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'
        assert output_handler.get_conformation_library_path(
            create_missing_folders=False) == \
            './DataLocal/Conformations/' \
            + '{}.conformation'.format(tag.upper()), \
            'Unexpected default conformation library path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, opls2005, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/OPLS2005/'
                             + 'HeteroAtoms/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'
            assert output_handler.get_conformation_library_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output/DataLocal/Conformations/' \
                + '{}.conformation'.format(tag.upper())), \
                'Unexpected default conformation library path'

    def test_datalocal_paths_for_offopls(self):
        """
        It tests the datalocal paths assignment for OpenFF-OPLS2005
        force field.
        """
        from peleffy.utils import OutputPathHandler
        from peleffy.forcefield import OpenFFOPLS2005ForceField

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')
        molecule._forcefield = OpenFFOPLS2005ForceField('OPLS2005')

        # Load force field
        hybridff = OpenFFOPLS2005ForceField(
            'openff_unconstrained-1.2.1.offxml')

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, hybridff,
                                           as_datalocal=True)

        # Validate output paths
        assert output_handler.get_rotamer_library_path(
            create_missing_folders=False) == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path(
            create_missing_folders=False) == \
            './DataLocal/Templates/OpenFF/Parsley/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path(
            create_missing_folders=False) == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'
        assert output_handler.get_conformation_library_path(
            create_missing_folders=False) == \
            './DataLocal/Conformations/' \
            + '{}.conformation'.format(tag.upper()), \
            'Unexpected default conformation library path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, hybridff, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/'
                             + 'OpenFF/Parsley/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'
            assert output_handler.get_conformation_library_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output/DataLocal/Conformations/' \
                + '{}.conformation'.format(tag.upper())), \
                'Unexpected default conformation library path'

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, hybridff,
                                           as_datalocal=True)

        # Set force field source for nonbonding parameters
        hybridff.set_nonbonding_parameters('opls2005')

        # Validate output paths
        assert output_handler.get_rotamer_library_path(
            create_missing_folders=False) == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path(
            create_missing_folders=False) == \
            './DataLocal/Templates/OPLS2005/HeteroAtoms/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path(
            create_missing_folders=False) == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'
        assert output_handler.get_conformation_library_path(
            create_missing_folders=False) == \
            './DataLocal/Conformations/' \
            + '{}.conformation'.format(tag.upper()), \
            'Unexpected default conformation library path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, hybridff, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/OPLS2005/'
                             + 'HeteroAtoms/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'
            assert output_handler.get_conformation_library_path(
                create_missing_folders=False) == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/Conformations/',
                             '{}.conformation'.format(tag.upper())), \
                'Unexpected default conformation library path'

    def test_folder_creation(self):
        """
        It tests the folder creation of the OutputPathHandler class.
        """
        from peleffy.utils import OutputPathHandler, temporary_cd
        from peleffy.forcefield import OpenForceField

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')

        # Load force field
        openff = OpenForceField('openff_unconstrained-1.2.1.offxml')

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            # Initialize output handler without output_path
            output_handler = OutputPathHandler(
                molecule, openff, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            # Test path getter without folder creation
            path = output_handler.get_rotamer_library_path(
                create_missing_folders=False)
            path_dir = os.path.dirname(path)
            assert os.path.isdir(path_dir) is False, \
                'This directory should not exist'
            path = output_handler.get_impact_template_path(
                create_missing_folders=False)
            path_dir = os.path.dirname(path)
            assert os.path.isdir(path_dir) is False, \
                'This directory should not exist'
            path = output_handler.get_solvent_template_path(
                create_missing_folders=False)
            path_dir = os.path.dirname(path)
            assert os.path.isdir(path_dir) is False, \
                'This directory should not exist'
            path = output_handler.get_conformation_library_path(
                create_missing_folders=False)
            path_dir = os.path.dirname(path)
            assert os.path.isdir(path_dir) is False, \
                'This directory should not exist'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                # Initialize output handler without output_path
                output_handler = OutputPathHandler(
                    molecule, openff, as_datalocal=True,
                    output_path=os.path.join(tmpdir, 'output'))

                # Test path getter with folder creation
                path = output_handler.get_rotamer_library_path(
                    create_missing_folders=True)
                path_dir = os.path.dirname(path)
                assert os.path.isdir(path_dir) is True, \
                    'This directory should exist'
                path = output_handler.get_impact_template_path(
                    create_missing_folders=True)
                path_dir = os.path.dirname(path)
                assert os.path.isdir(path_dir) is True, \
                    'This directory should exist'
                path = output_handler.get_solvent_template_path(
                    create_missing_folders=True)
                path_dir = os.path.dirname(path)
                assert os.path.isdir(path_dir) is True, \
                    'This directory should exist'
                path = output_handler.get_conformation_library_path(
                    create_missing_folders=True)
                path_dir = os.path.dirname(path)
                assert os.path.isdir(path_dir) is True, \
                    'This directory should exist'


class TestMAEParser(object):
    """
    It contains all the tests to validate functions that parse MAE files.
    """

    def test_parse_charges_from_mae(self):
        """
        It tests the MAE parser for assigning partial charges from an external
        file.
        """
        from peleffy.utils import get_data_file_path, parse_charges_from_mae
        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenForceField

        FORCEFIELD = 'openff_unconstrained-1.2.0.offxml'

        PATH_MAE_MAT = get_data_file_path('tests/MAT.mae')
        PATH_PDB_MAT = get_data_file_path('tests/MAT.pdb')
        PATH_MAE_BHP = get_data_file_path('ligands/BHP.mae')
        PATH_PDB_BHP = get_data_file_path('ligands/BHP.pdb')
        PATH_MAE_ETL = get_data_file_path('ligands/ethylene.mae')
        PATH_PDB_ETL = get_data_file_path('ligands/ethylene.pdb')
        PATH_PDB_MAL = get_data_file_path('ligands/malonate.pdb')

        CHARGES_REFERENCE_BHP = [[-0.35703, 'O1'],
                                 [-0.59535, 'O2'],
                                 [-0.50292, 'O3'],
                                 [-0.25243, 'C1'],
                                 [0.30438, 'C2'],
                                 [0.22092, 'C3'],
                                 [0.54336, 'C4'],
                                 [-0.20569, 'C5'],
                                 [-0.20192, 'C6'],
                                 [0.16631, 'C7'],
                                 [0.02422, 'C8'],
                                 [-0.09115, 'C9'],
                                 [-0.09904, 'C10'],
                                 [-0.15673, 'C11'],
                                 [-0.13245, 'C12'],
                                 [-0.17806, 'C13'],
                                 [-0.12489, 'C14'],
                                 [-0.09307, 'C15'],
                                 [-0.08973, 'C16'],
                                 [0.05397, 'H1'],
                                 [0.07338, 'H2'],
                                 [0.04514, 'H3'],
                                 [0.12979, 'H4'],
                                 [0.11025, 'H5'],
                                 [0.04054, 'H6'],
                                 [0.04581, 'H7'],
                                 [0.11444, 'H8'],
                                 [0.11761, 'H9'],
                                 [0.37274, 'H10'],
                                 [0.11825, 'H11'],
                                 [0.12584, 'H12'],
                                 [0.1381, 'H13'],
                                 [0.11696, 'H14'],
                                 [0.11148, 'H15'],
                                 [0.10697, 'H16']]

        CHARGES_REFERENCE_MAT = [[-0.18938, 'F1'],
                                 [-0.21715, 'F2'],
                                 [-0.21234, 'F3'],
                                 [-0.39736, 'O1'],
                                 [-0.58890, 'O2'],
                                 [-1.00825, 'N1'],
                                 [0.72066, 'C1'],
                                 [-0.06281, 'C2'],
                                 [-0.67474, 'C3'],
                                 [0.10391, 'C4'],
                                 [0.16293, 'C5'],
                                 [-0.61076, 'C6'],
                                 [0.78183, 'C7'],
                                 [0.27041, 'C8'],
                                 [-0.48769, 'C9'],
                                 [0.15704, 'C10'],
                                 [-0.02646, 'H1'],
                                 [-0.08394, 'H2'],
                                 [-0.01308, 'H3'],
                                 [0.14006, 'H4'],
                                 [0.12960, 'H5'],
                                 [0.31245, 'H6'],
                                 [0.30268, 'H7'],
                                 [0.17026, 'H8'],
                                 [0.15782, 'H9'],
                                 [0.03175, 'C11'],
                                 [0.03894, 'H10'],
                                 [0.07509, 'H11'],
                                 [0.01743, 'H12']]

        # Check up correct charges for malonate
        m = Molecule(PATH_PDB_MAT)
        openff = OpenForceField(FORCEFIELD)
        parameters = openff.parameterize(m, charge_method='dummy')
        parameters = parse_charges_from_mae(PATH_MAE_MAT, parameters)

        for charge, atom_name in zip(parameters['charges'],
                                     parameters['atom_names']):
            assert [charge._value, atom_name.replace('_', '')] in \
                CHARGES_REFERENCE_MAT, \
                'Incorrect charge value for {}.'.format(atom_name)

        # Check up correct charges for BHP
        m = Molecule(PATH_PDB_BHP)
        openff = OpenForceField(FORCEFIELD)
        parameters = openff.parameterize(m, charge_method='dummy')
        parameters = parse_charges_from_mae(PATH_MAE_BHP, parameters)

        for charge, atom_name in zip(parameters['charges'],
                                     parameters['atom_names']):
            assert [charge._value, atom_name.replace('_', '')] in \
                CHARGES_REFERENCE_BHP, \
                'Incorrect charge value for {}.'.format(atom_name)

        # Error: MAE file without charges information
        m = Molecule(PATH_PDB_ETL)
        openff = OpenForceField(FORCEFIELD)
        parameters = openff.parameterize(m, charge_method='dummy')
        with pytest.raises(ValueError):
            _ = parse_charges_from_mae(PATH_MAE_ETL, parameters)

        # Error: Inconsistency between Moelcule atom names and MAE atom names
        m = Molecule(PATH_PDB_MAL)
        openff = OpenForceField(FORCEFIELD)
        parameters = openff.parameterize(m, charge_method='dummy')
        with pytest.raises(ValueError):
            _ = parse_charges_from_mae(PATH_PDB_BHP, parameters)


class TestPDBFile(object):
    """
    It contains all the tests to validate the PDBFile class.
    """
    def test_get_molecule_from_chain(self):
        """
        It tests the method that returns a molecule from a selected chain from
        an input PDB file.
        """
        def compare_molecules(mol1, mol2):
            """
            It compares that two molecules contain the same PDB atom names and
            that have been initialized with the same flags.

            Parameters
            ----------
            mol1 : a peleffy.topology.Molecule object
                Molecule.
            mol2 : a peleffy.topology.Molecule object
                Reference molecule.

            Returns
            -------
            are_equal : bool
                Wheter the molecules are equal or not.
            """
            if mol1._rotamer_resolution == mol2._rotamer_resolution \
                    and mol1._exclude_terminal_rotamers ==  \
            mol2._exclude_terminal_rotamers\
                    and mol1._allow_undefined_stereo == \
            mol2._allow_undefined_stereo\
                    and mol1.get_pdb_atom_names() == \
            mol2.get_pdb_atom_names():
                return True
            else:
                return False

        from peleffy.utils.input import PDBFile
        from peleffy.utils import get_data_file_path
        from peleffy.topology import Molecule

        PATH_COMPLEX_PDB = get_data_file_path('complexes/LYS_BNZ.pdb')
        PATH_LIGAND_PDB = get_data_file_path('ligands/BNZ.pdb')

        # Test method get_molecule_from_chain
        PDBreader = PDBFile(PATH_COMPLEX_PDB)
        molecules = PDBreader.get_molecules_from_chain(selected_chain='L')
        assert len(molecules) == 1, 'Unexpected number of molecules'
        assert compare_molecules(molecules[0], Molecule(PATH_LIGAND_PDB)) is True

        # Test allow_undefined_stereo flag
        PDBreader = PDBFile(PATH_COMPLEX_PDB)
        molecules = \
            PDBreader.get_molecules_from_chain(selected_chain='L',
                                               allow_undefined_stereo=True)
        assert len(molecules) == 1, 'Unexpected number of molecules'
        assert compare_molecules(molecules[0],
                                 Molecule(PATH_LIGAND_PDB,
                                          allow_undefined_stereo=True)) is True

        # Test exclude_terminal_rotamers flag
        PDBreader = PDBFile(PATH_COMPLEX_PDB)
        molecules = \
            PDBreader.get_molecules_from_chain(selected_chain='L',
                                               exclude_terminal_rotamers=False)
        assert len(molecules) == 1, 'Unexpected number of molecules'
        assert compare_molecules(molecules[0], Molecule(
            PATH_LIGAND_PDB, exclude_terminal_rotamers=False)) is True

        # Test rotamer_resolution flag
        PDBreader = PDBFile(PATH_COMPLEX_PDB)
        molecules = PDBreader.get_molecules_from_chain(selected_chain='L',
                                                       rotamer_resolution=10)
        assert compare_molecules(molecules[0],
                                 Molecule(PATH_LIGAND_PDB,
                                          rotamer_resolution=10)) is True

        # Test core_constraints flag
        PDBreader = PDBFile(PATH_COMPLEX_PDB)
        molecules = PDBreader.get_molecules_from_chain(selected_chain='L',
                                                       core_constraints=[1,2])
        assert len(molecules) == 1, 'Unexpected number of molecules'
        assert molecules[0].core_constraints == [1,2]


    def test_get_hetero_molecules(self):
        """
        It tests the method that returns a list with all the hetero molecules in
        an input PDB.
        """
        def compare_molecules(mol1, mol2):
            """
            It compares that two molecules contain the same PDB atom names and
            that have been initialized with the same flags.

            Parameters
            ----------
            mol1 : a peleffy.topology.Molecule object
                Molecule.
            mol2 : a peleffy.topology.Molecule object
                Reference molecule.

            Returns
            -------
            are_equal : bool
                Wheter the molecules are equal or not.
            """
            if mol1._rotamer_resolution == mol2._rotamer_resolution \
                    and mol1._exclude_terminal_rotamers ==  \
            mol2._exclude_terminal_rotamers\
                    and mol1._allow_undefined_stereo == \
            mol2._allow_undefined_stereo\
                    and mol1.get_pdb_atom_names() == \
            mol2.get_pdb_atom_names():
                return True
            else:
                return False

        from peleffy.utils.input import PDBFile
        from peleffy.utils import get_data_file_path
        from peleffy.topology import Molecule

        PATH_COMPLEX_PDB = get_data_file_path('complexes/LIG1_LIG2.pdb')
        PATH_LIGAND1 = get_data_file_path('ligands/LIG1.pdb')
        PATH_LIGAND2 = get_data_file_path('ligands/LIG2.pdb')

        # Test extract multiple molecules
        PDBreader = PDBFile(PATH_COMPLEX_PDB)
        molecules = PDBreader.get_hetero_molecules()

        for molecule in molecules:
            for molecule_ref in [Molecule(PATH_LIGAND1),Molecule(PATH_LIGAND2)]:
                if molecule._tag == molecule_ref.tag:
                    assert compare_molecules(molecule, molecule_ref)

        # Test set ligand core constraints
        PATH_COMPLEX_PDB = get_data_file_path('complexes/complex_test.pdb')
        molecules = \
            PDBreader.get_hetero_molecules(allow_undefined_stereo = True,
                                           ligand_core_constraints = [1,2,3],
                                           ligand_resname = 'ANL')
        for molecule in molecules:
            if molecule._tag == 'ANL':
                assert molecule.core_constraints == [1,2,3]
            else:
                assert molecule.core_constraints == []


    def test_raise_errors(self):
        """
        It tests all the possible errors when geting a molecule from a chain.
        """
        from peleffy.utils.input import PDBFile
        from peleffy.utils import get_data_file_path

        PATH_COMPLEX_PDB = get_data_file_path('complexes/LYS_BNZ.pdb')

        # The chain selected does not exist in the PDB file
        with pytest.raises(ValueError):
            PDBreader = PDBFile(PATH_COMPLEX_PDB)
            _ = PDBreader.get_molecules_from_chain(selected_chain='F')

        # The chain selected is not an hetero molecule
        with pytest.raises(ValueError):
            PDBreader = PDBFile(PATH_COMPLEX_PDB)
            _ = PDBreader.get_molecules_from_chain(selected_chain='A')
