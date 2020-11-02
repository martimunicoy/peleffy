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

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')
        molecule.parameterize('openff_unconstrained-1.2.1.offxml',
                              charge_method='gasteiger')

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, as_datalocal=False)

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

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, as_datalocal=False,
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

    def test_datalocal_paths_for_openff(self):
        """It tests the datalocal paths assignment for OpenFF."""
        from peleffy.utils import OutputPathHandler

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')
        molecule.parameterize('openff_unconstrained-1.2.1.offxml',
                              charge_method='gasteiger')

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, as_datalocal=True)

        # Validate output paths
        assert output_handler.get_rotamer_library_path() == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path() == \
            './DataLocal/Templates/OFF/Parsley/HeteroAtoms/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path() == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/OFF/Pars'
                             + 'ley/HeteroAtoms/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path() == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'

    def test_datalocal_paths_for_opls(self):
        """It tests the datalocal paths assignment for OPLS2005."""
        from peleffy.utils import OutputPathHandler
        from peleffy.forcefield import OPLS2005ForceField

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')
        molecule._forcefield = OPLS2005ForceField('OPLS2005')

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, as_datalocal=True)

        # Validate output paths
        assert output_handler.get_rotamer_library_path() == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path() == \
            './DataLocal/Templates/OPLS2005/HeteroAtoms/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path() == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/OPLS2005/'
                             + 'HeteroAtoms/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path() == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'

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

        # Molecule's tag
        tag = molecule.tag

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, as_datalocal=True)

        # Validate output paths
        assert output_handler.get_rotamer_library_path() == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path() == \
            './DataLocal/Templates/OFF/Parsley/HeteroAtoms/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path() == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/OFF/Pars'
                             + 'ley/HeteroAtoms/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path() == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'

        # Initialize output handler without output_path
        output_handler = OutputPathHandler(molecule, as_datalocal=True)

        # Set force field source for nonbonding parameters
        molecule._forcefield.set_nonbonding_parameters('opls2005')

        # Validate output paths
        assert output_handler.get_rotamer_library_path() == \
            './DataLocal/LigandRotamerLibs/' \
            + '{}.rot.assign'.format(tag.upper()), \
            'Unexpected default rotamer library path'
        assert output_handler.get_impact_template_path() == \
            './DataLocal/Templates/OPLS2005/HeteroAtoms/' \
            + '{}z'.format(tag.lower()), \
            'Unexpected default Impact template path'
        assert output_handler.get_solvent_template_path() == \
            './DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            output_handler = OutputPathHandler(
                molecule, as_datalocal=True,
                output_path=os.path.join(tmpdir, 'output'))

            assert output_handler.get_rotamer_library_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/LigandRotamerLibs/'
                             + '{}.rot.assign'.format(tag.upper())), \
                'Unexpected default rotamer library path'
            assert output_handler.get_impact_template_path() == \
                os.path.join(tmpdir, 'output', 'DataLocal/Templates/OPLS2005/'
                             + 'HeteroAtoms/{}z'.format(tag.lower())), \
                'Unexpected default Impact template path'
            assert output_handler.get_solvent_template_path() == \
                os.path.join(tmpdir, 'output',
                             'DataLocal/OBC/ligandParams.txt'), \
                'Unexpected default solvent parameters path'

    def test_folder_creation(self):
        """
        It tests the folder creation of the OutputPathHandler class.
        """
        from peleffy.utils import OutputPathHandler

        # Load benzene molecule
        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')
        molecule.parameterize('openff_unconstrained-1.2.1.offxml',
                              charge_method='gasteiger')

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            # Initialize output handler without output_path
            output_handler = OutputPathHandler(
                molecule, as_datalocal=True,
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

        # Initialize output handler with an output_path set
        with tempfile.TemporaryDirectory() as tmpdir:
            # Initialize output handler without output_path
            output_handler = OutputPathHandler(
                molecule, as_datalocal=True,
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
