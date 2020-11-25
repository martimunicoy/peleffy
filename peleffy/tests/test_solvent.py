"""
This module contains the tests to check peleffy's solvent.
"""

import tempfile

from peleffy.utils import get_data_file_path, temporary_cd
from peleffy.topology import Molecule, Topology
from peleffy.forcefield import OPLS2005ForceField
from peleffy.solvent import OPLSOBC


class TestSolvent(object):
    """
    It contains all the tests that validate the solvent-template generator.
    """

    def test_OBCOPLS_writer(self):
        """
        It test the function that writes a OPLS2005CompatibleSolvent object to
        a file compatible with PELE.
        """
        from .utils import parameterize_opls2005, compare_files_without_order

        TEMPLATE_PARAMS_ETL = get_data_file_path(
            'tests/ETL_solventParamsHCTOBC.txt')
        TEMPLATE_PARAMS_MAL = get_data_file_path(
            'tests/MAL_solventParamsHCTOBC.txt')
        TEMPLATE_PARAMS_MET = get_data_file_path(
            'tests/MET_solventParamsHCTOBC.txt')

        def test_OBCOPLS_writer_ligand(pdbfile, tag_name, ffld_name,
                                       reference_file):
            """
            Given a ligand, it tests that the output parameters file
            corresponds to the refenrece file.

            Parameters
            ----------
            pdbfile : str
                The path to the PDB of the ligand to test
            ffld_name : str
                The path to the ffld_server's output file
            reference_file : str
                The path to reference TXT file compatible with PELE
            """
            with tempfile.TemporaryDirectory() as tmpdir:
                with temporary_cd(tmpdir):

                    # Loads the  molecule
                    molecule = Molecule(get_data_file_path(pdbfile),
                                        tag=tag_name)

                    # Sets forcefield and parameterizes it
                    opls2005 = OPLS2005ForceField()
                    ffld_file = get_data_file_path(ffld_name)
                    parameters = parameterize_opls2005(opls2005,
                                                       molecule,
                                                       ffld_file)

                    # Initializes topology
                    topology = Topology(molecule, parameters)

                    # Initializes solvent and getsparameters file
                    solvent = OPLSOBC(topology)
                    solvent.to_file('OBC_parameters.txt')

                    # Compare the output file with the reference parameters file
                    compare_files_without_order('OBC_parameters.txt',
                                                reference_file)

        # Test for ethylene
        test_OBCOPLS_writer_ligand(pdbfile='ligands/ethylene.pdb',
                                   tag_name='ETL',
                                   ffld_name='tests/ETL_ffld_output.txt',
                                   reference_file=TEMPLATE_PARAMS_ETL)

        # Test for methane
        test_OBCOPLS_writer_ligand(pdbfile='ligands/methane.pdb',
                                   tag_name='MET',
                                   ffld_name='tests/MET_ffld_output.txt',
                                   reference_file=TEMPLATE_PARAMS_MET)

        # Test for malonate
        test_OBCOPLS_writer_ligand(pdbfile='ligands/malonate.pdb',
                                   tag_name='MAL',
                                   ffld_name='tests/MAL_ffld_output.txt',
                                   reference_file=TEMPLATE_PARAMS_MAL)
