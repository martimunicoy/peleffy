"""
This module contains the tests to check offpele's molecular representations.
"""

import pytest

import os
import tempfile
from offpele.main import run_offpele, handle_output_paths
from offpele.utils import get_data_file_path, temporary_cd
from offpele.topology import Molecule


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'


class TestMain(object):
    """
    It wraps all tests that involve the Molecule class.
    """

    def test_offpele_default_call(self):
        """
        It checks the default call of offpele's main function.
        """
        LIGAND_PATH = 'ligands/BNZ.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                run_offpele(ligand_path, output=tmpdir)

    def test_offpele_custom_call(self):
        """
        It checks the custom call of offpele's main function.
        """
        LIGAND_PATH = 'ligands/BNZ.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                run_offpele(ligand_path,
                            forcefield=FORCEFIELD_NAME,
                            resolution=10,
                            charges_method='gasteiger',
                            output=tmpdir,
                            with_solvent=True,
                            as_datalocal=True)

    def test_default_output_paths(self):
        """
        It checks the default output paths that are used for each parameter
        file from offpele.
        """

        def from_PosixPath_to_string(paths):
            """
            Convert PosixPaths to strings
            """
            return map(str, paths)

        molecule = Molecule(smiles='c1ccccc1', name='benzene', tag='BNZ')

        rotlib_path, impact_path, solvent_path = \
            handle_output_paths(molecule, '', False)

        # Convert PosixPaths to strings
        rotlib_path, impact_path, solvent_path = map(
            str, [rotlib_path, impact_path, solvent_path])

        assert rotlib_path == 'BNZ.rot.assign', 'Unexpected default ' \
            + 'rotamer library path'
        assert impact_path == 'bnzz', 'Unexpected default Impact ' \
            + 'template path'
        assert solvent_path == 'ligandParams.txt', 'Unexpected default ' \
            + 'solvent parameters path'

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                # To avoid the complain about unexistent folder
                os.mkdir('output')
                rotlib_path, impact_path, solvent_path = \
                    handle_output_paths(molecule, 'output', False)

        # Convert PosixPaths to strings
        rotlib_path, impact_path, solvent_path = map(
            str, [rotlib_path, impact_path, solvent_path])

        assert rotlib_path == 'output/BNZ.rot.assign', 'Unexpected default ' \
            + 'rotamer library path'
        assert impact_path == 'output/bnzz', 'Unexpected default Impact ' \
            + 'template path'
        assert solvent_path == 'output/ligandParams.txt', 'Unexpected ' \
            + 'default solvent parameters path'

        rotlib_path, impact_path, solvent_path = \
            handle_output_paths(molecule, '', True)

        # Convert PosixPaths to strings
        rotlib_path, impact_path, solvent_path = map(
            str, [rotlib_path, impact_path, solvent_path])

        assert rotlib_path == 'DataLocal/LigandRotamerLibs/' \
            + 'BNZ.rot.assign', 'Unexpected default rotamer library path'
        assert impact_path == 'DataLocal/Templates/OFF/Parsley/' \
            + 'HeteroAtoms/bnzz', 'Unexpected default Impact template'
        assert solvent_path == 'DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                # To avoid the complain about unexistent folder
                os.mkdir('output')
                rotlib_path, impact_path, solvent_path = \
                    handle_output_paths(molecule, 'output', True)

        # Convert PosixPaths to strings
        rotlib_path, impact_path, solvent_path = map(
            str, [rotlib_path, impact_path, solvent_path])

        assert rotlib_path == 'output/DataLocal/LigandRotamerLibs/' \
            + 'BNZ.rot.assign', 'Unexpected default rotamer library path'
        assert impact_path == 'output/DataLocal/Templates/OFF/Parsley/' \
            + 'HeteroAtoms/bnzz', 'Unexpected default Impact template path'
        assert solvent_path == 'output/DataLocal/OBC/ligandParams.txt', \
            'Unexpected default solvent parameters path'
