"""
This module contains the tests to check peleffy's solvent.
"""

import tempfile

from peleffy.utils import get_data_file_path, temporary_cd
from peleffy.topology import Molecule, Topology
from peleffy.forcefield import OPLS2005ForceField, OpenForceField
from peleffy.solvent import OPLSOBC, OBC2


class TestSolvent(object):
    """
    It contains all the tests that validate the solvent-template generator.
    """

    def test_single_topology(self):
        """
        It tests the class that generates a OpenFFCompatibleSolvent object for
        a single topology.
        """
        from .utils import compare_dicts
        import json

        TEMPLATE_PARAMS_MAL = get_data_file_path('tests/ligandParams_MAL.txt')

        # Loads the  molecule
        molecule = Molecule(path=get_data_file_path('ligands/malonate.pdb'),
                            tag='MAL')

        # Sets forcefield and parameterizes it
        ff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Initializes topology
        topology = Topology(molecule, parameters)

        # Initializes solvent and gets parameters file
        solvent = OBC2(topology)
        solvent_dict = solvent.to_dict()

        # Loads reference dict from template
        with open(TEMPLATE_PARAMS_MAL, 'r') as f:
            reference_dict = json.load(f)

        # Compare the output parameters dict with the reference parameters
        compare_dicts(reference_dict, solvent_dict)

    def test_multiple_topologies(self):
        """
        It tests the class that generates a OpenFFCompatibleSolvent object for
        multiple topologies.
        """
        from .utils import compare_dicts, merge_dicts

        # Path to multiple non standard residues
        pdb_path_MAL = get_data_file_path('ligands/malonate.pdb')
        pdb_path_MET = get_data_file_path('ligands/methane.pdb')

        # Force Field to parameterize the molecules
        ff = OpenForceField('openff_unconstrained-1.2.1.offxml')

        # Topology of malonate
        mol_MAL = Molecule(path=pdb_path_MAL, tag='MAL')
        parameters_MAL = ff.parameterize(mol_MAL, charge_method='gasteiger')
        topology_MAL = Topology(mol_MAL, parameters_MAL)

        # Topology of methane
        mol_MET = Molecule(path=pdb_path_MET, tag='MET')
        parameters_MET = ff.parameterize(mol_MET, charge_method='gasteiger')
        topology_MET = Topology(mol_MET, parameters_MET)

        # List containing both topologies
        topologies = [topology_MAL, topology_MET]

        # Generate the Solvent parameters dictionaries
        solvent_MAL_dict = OBC2(topology_MAL).to_dict()
        solvent_MET_dict = OBC2(topology_MET).to_dict()
        solvent_dict = OBC2(topologies).to_dict()

        # Check that merging both single topology dicitionaries we obtain the
        # same dictionary that using multiple topologies
        compare_dicts(merge_dicts(solvent_MAL_dict['SolventParameters'],
                                  solvent_MET_dict['SolventParameters']),
                      solvent_dict['SolventParameters'])

    def test_multiple_topologies_writer(self):
        """
        It tests the class that generates a OpenFFCompatibleSolvent object for multiple topologies. It compares the outcome of the Solvent writer with
        a reference file.
        """
        from .utils import compare_dicts, parameterize_opls2005
        import json

        TEMPLATE_PARAMS = get_data_file_path('tests/ligandParams.txt')

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                path_OXO = get_data_file_path('tests/MRO_oleic/OXO.pdb')
                path_OLC = get_data_file_path('tests/MRO_oleic/OLC.pdb')

                ff = OpenForceField('openff_unconstrained-1.2.1.offxml')
                opls2005 = OPLS2005ForceField()

                # Group OXO
                m_OXO = Molecule(path_OXO)
                ffld_file = get_data_file_path('tests/OXO_ffld_output.txt')
                parameters_OXO = parameterize_opls2005(opls2005, m_OXO,
                                                       ffld_file)
                topology_OXO = Topology(m_OXO, parameters_OXO)

                # Acid oleic
                m_OLC = Molecule(path_OLC)
                parameters_OLC = ff.parameterize(m_OLC,
                                                 charge_method='gasteiger')
                topology_OLC = Topology(m_OLC, parameters_OLC)

                # Multiple topologies
                topologies = [topology_OXO, topology_OLC]
                solvent = OBC2(topologies)
                solvent.to_file('OBC_parameters.txt')

                # Loads reference dict from template
                with open(TEMPLATE_PARAMS, 'r') as f:
                    reference_dict = json.load(f)

                # Loads the generated template into a dict
                with open('OBC_parameters.txt', 'r') as f:
                    solvent_dict = json.load(f)

                # Compare the output parameters dict with the reference parameters
                compare_dicts(reference_dict, solvent_dict)

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

                    # Initializes solvent and gets parameters file
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
