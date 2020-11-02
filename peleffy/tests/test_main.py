"""
This module contains the tests to check peleffy's molecular representations.
"""

import pytest

import tempfile
from peleffy.utils import get_data_file_path, temporary_cd


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'


class TestMain(object):
    """
    It wraps all tests that involve the Molecule class.
    """

    def test_peleffy_default_call(self):
        """
        It checks the default call of peleffy's main function.
        """
        from peleffy.main import run_peleffy

        LIGAND_PATH = 'ligands/BNZ.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                run_peleffy(ligand_path, output=tmpdir)

    def test_peleffy_custom_call(self):
        """
        It checks the custom call of peleffy's main function.
        """
        from peleffy.main import run_peleffy

        LIGAND_PATH = 'ligands/BNZ.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                run_peleffy(ligand_path,
                            forcefield=FORCEFIELD_NAME,
                            resolution=10,
                            charge_method='gasteiger',
                            output=tmpdir,
                            with_solvent=True,
                            as_datalocal=True)

    def test_peleffy_argparse(self):
        """It checks the command-line argument parser of peleffy."""
        from peleffy.main import parse_args

        # Test positional arguments
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            parsed_args = parse_args([])

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

        # Test defaults
        parsed_args = parse_args(['BNZ.pdb'])

        assert parsed_args.as_datalocal is False, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'am1bcc', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is False, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.2.0.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is False, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output is None, \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'BNZ.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 30, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is False, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        # Test custom shorts
        parsed_args = parse_args(['TOL.pdb',
                                  '-f', 'openff_unconstrained-1.0.0.offxml',
                                  '-r', '60',
                                  '-o', 'my_custom_output',
                                  '-c', 'gasteiger'])

        assert parsed_args.as_datalocal is False, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'gasteiger', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is False, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.0.0.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is False, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output == 'my_custom_output', \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'TOL.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 60, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is False, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        # Test custom longs
        parsed_args = parse_args(['MET.pdb',
                                  '--forcefield',
                                  'openff_unconstrained-1.0.1.offxml',
                                  '--resolution', '120',
                                  '--output', 'my_custom_output2',
                                  '--charge_method', 'gasteiger'])

        assert parsed_args.as_datalocal is False, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'gasteiger', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is False, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.0.1.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is False, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output == 'my_custom_output2', \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'MET.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 120, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is False, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        # Test unexpected charge method
        with pytest.raises(SystemExit) as pytest_wrapped_e:
            parsed_args = parse_args(['TOL.pdb', '-c', 'unexpected'])

        assert pytest_wrapped_e.type == SystemExit
        assert pytest_wrapped_e.value.code == 2

        # Test as_datalocal argument
        parsed_args = parse_args(['MET.pdb',
                                  '--as_datalocal'])

        assert parsed_args.as_datalocal is True, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'am1bcc', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is False, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.2.0.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is False, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output is None, \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'MET.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 30, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is False, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        # Test include_terminal_rotamers argument
        parsed_args = parse_args(['MET.pdb',
                                  '--include_terminal_rotamers'])

        assert parsed_args.as_datalocal is False, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'am1bcc', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is False, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.2.0.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is True, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output is None, \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'MET.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 30, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is False, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        # Test silent argument
        parsed_args = parse_args(['MET.pdb',
                                  '-s'])

        assert parsed_args.as_datalocal is False, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'am1bcc', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is False, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.2.0.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is False, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output is None, \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'MET.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 30, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is True, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        parse_args(['MET.pdb', '-s']) == parse_args(['MET.pdb', '--silent'])

        # Test debug argument
        parsed_args = parse_args(['MET.pdb',
                                  '-d'])

        assert parsed_args.as_datalocal is False, \
            'Unexpected as_datalocal settings were parsed'
        assert parsed_args.charge_method == 'am1bcc', \
            'Unexpected charge_method settings were parsed'
        assert parsed_args.debug is True, \
            'Unexpected debug settings were parsed'
        assert parsed_args.forcefield == 'openff_unconstrained-1.2.0.offxml', \
            'Unexpected forcefield settings were parsed'
        assert parsed_args.include_terminal_rotamers is False, \
            'Unexpected include_terminal_rotamers settings were parsed'
        assert parsed_args.output is None, \
            'Unexpected output settings were parsed'
        assert parsed_args.pdb_file == 'MET.pdb', \
            'Unexpected pdb_file settings were parsed'
        assert parsed_args.resolution == 30, \
            'Unexpected resolution settings were parsed'
        assert parsed_args.silent is False, \
            'Unexpected silent settings were parsed'
        assert parsed_args.with_solvent is False, \
            'Unexpected with_solvent settings were parsed'

        parse_args(['MET.pdb', '-d']) == parse_args(['MET.pdb', '--debug'])

    def test_peleffy_main(self):
        """It checks the main function of peleffy."""
        from peleffy.main import parse_args, main
        from peleffy.utils import Logger
        import logging

        ligand_path = get_data_file_path('ligands/BNZ.pdb')

        # Test default settings
        args = parse_args([ligand_path])
        main(args)

        logger = Logger()
        for handler in logger._logger.handlers:
            assert handler.level == logging.INFO

        # Test silent settings
        args = parse_args([ligand_path, '--silent'])
        main(args)

        logger = Logger()
        for handler in logger._logger.handlers:
            assert handler.level == logging.CRITICAL

        # Test silent settings
        args = parse_args([ligand_path, '--debug'])
        main(args)

        logger = Logger()
        for handler in logger._logger.handlers:
            assert handler.level == logging.DEBUG
