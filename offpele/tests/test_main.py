"""
This module contains the tests to check offpele's molecular representations.
"""

import pytest

import tempfile
from offpele.main import run_offpele
from offpele.utils import get_data_file_path, temporary_cd


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'


class TestMain(object):
    """
    It wraps all tests that involve the Molecule class.
    """

    def test_offpele_default_call(self):
        LIGAND_PATH = 'ligands/BNZ.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                run_offpele(ligand_path, output=tmpdir)

    def test_offpele_custom_call(self):
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
