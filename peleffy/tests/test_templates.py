"""
This module contains the tests to check the peleffy's Impact template builder.
"""

import pytest

import tempfile

from peleffy.utils import get_data_file_path, temporary_cd
from peleffy.topology import Molecule
from peleffy.template import Impact
from peleffy.tests.utils import compare_files
from peleffy.topology import Topology
from peleffy.forcefield import OpenForceField, OPLS2005ForceField


class TestImpactTemplate(object):
    """It tests the Impact class that writes a PELE's Impact template"""

    OPENFF_FORCEFIELD = 'openff_unconstrained-1.2.1.offxml'

    def test_input(self):
        """
        It tests that the topology given to Impact() is of the correct
        format, peleffy.topology.Topology.
        """
        from peleffy.forcefield.parameters import BaseParameterWrapper

        LIGAND_PATH = 'ligands/benzene.pdb'

        ligand_path = get_data_file_path(LIGAND_PATH)
        molecule = Molecule(ligand_path)

        parameters = BaseParameterWrapper()

        topology = Topology(molecule, parameters)

        # Impact() gets nothing as argument
        with pytest.raises(TypeError):
            _ = Impact()

        # Impact() gets a non Topology object
        with pytest.raises(TypeError):
            _ = Impact('passing a str instead of a Topology')

        # This should work
        _ = Impact(topology)

    def test_writer_OFF(self):
        """
        It tests the writer attribute of the Impact class using OFF
        to parameterize.
        """

        TEMPLATE_METZ = get_data_file_path('tests/metz')
        TEMPLATE_MATZ = get_data_file_path('tests/malz')
        TEMPLATE_ETLZ = get_data_file_path('tests/etlz')

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                # Generates the template for methane
                pdb_path = get_data_file_path('ligands/methane.pdb')
                molecule = Molecule(pdb_path)
                openff = OpenForceField(self.OPENFF_FORCEFIELD)
                parameters = openff.parameterize(molecule)
                topology = Topology(molecule, parameters)

                # Generates the impact template for methane
                impact = Impact(topology)
                impact.to_file('metz')

                # Compare the reference template and the generated template
                compare_files(file1=TEMPLATE_METZ, file2='metz')

                # Generates the template for malonate
                pdb_path = get_data_file_path('ligands/malonate.pdb')
                molecule = Molecule(pdb_path)
                openff = OpenForceField(self.OPENFF_FORCEFIELD)
                parameters = openff.parameterize(molecule)
                topology = Topology(molecule, parameters)

                # Generates the impact template for malonate
                impact = Impact(topology)
                impact.to_file('malz')

                # Compare the reference template and the generated template
                compare_files(file1=TEMPLATE_MATZ, file2='malz')

                # Generates the template for ethylene
                pdb_path = get_data_file_path('ligands/ethylene.pdb')
                molecule = Molecule(pdb_path, tag='ETL')  # Note that in this case we are assigning a tag to the molecule which will be used in the Impact template
                openff = OpenForceField(self.OPENFF_FORCEFIELD)
                parameters = openff.parameterize(molecule)
                topology = Topology(molecule, parameters)

                # Generates the impact template for ethylene
                impact = Impact(topology)
                impact.to_file('etlz')

                # Compare the reference template and the generated template
                compare_files(file1=TEMPLATE_ETLZ, file2='etlz')

    def test_writer_OPLS(self):
        """
        It tests the writer attribute of the Impact class using OPLS to parameterize.
        """
        from .utils import parameterize_opls2005

        TEMPLATE_METZ_OPLS = get_data_file_path('tests/OPLS_metz')
        TEMPLATE_MALZ_OPLS = get_data_file_path('tests/OPLS_malz')
        TEMPLATE_ETLZ_OPLS = get_data_file_path('tests/OPLS_etlz')

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                # Generates the template for methane using OPLS
                opls2005 = OPLS2005ForceField()
                pdb_path = get_data_file_path('ligands/methane.pdb')
                molecule = Molecule(pdb_path)
                ffld_file = get_data_file_path('tests/MET_ffld_output.txt')
                parameters = parameterize_opls2005(opls2005,
                                                   molecule,
                                                   ffld_file)
                topology = Topology(molecule, parameters)

                # Generates the impact template for methane
                impact = Impact(topology)
                impact.to_file('metz')

                # Compare the reference template and the generated template
                compare_files(file1=TEMPLATE_METZ_OPLS, file2='metz')

                # Generates the template for malonate using OPLS
                opls2005 = OPLS2005ForceField()
                pdb_path = get_data_file_path('ligands/malonate.pdb')
                molecule = Molecule(pdb_path)
                ffld_file = get_data_file_path('tests/MAL_ffld_output.txt')
                parameters = parameterize_opls2005(opls2005,
                                                   molecule,
                                                   ffld_file)
                topology = Topology(molecule, parameters)

                # Generates the impact template for malonate
                impact = Impact(topology)
                impact.to_file('malz')

                # Compare the reference template and the generated template
                compare_files(file1=TEMPLATE_MALZ_OPLS, file2='malz')

                # Generates the template for ethylene using OPLS
                opls2005 = OPLS2005ForceField()
                pdb_path = get_data_file_path('ligands/ethylene.pdb')
                molecule = Molecule(pdb_path, tag='ETL')
                ffld_file = get_data_file_path('tests/ETL_ffld_output.txt')
                parameters = parameterize_opls2005(opls2005,
                                                   molecule,
                                                   ffld_file)
                topology = Topology(molecule, parameters)

                # Generates the impact template for ethylene
                impact = Impact(topology)
                impact.to_file('etlz')

                # Compare the reference template and the generated template
                compare_files(file1=TEMPLATE_ETLZ_OPLS, file2='etlz')

    def test_get_absolute_parent_atom(self):
        """
        It tests the _get_absolute_parent_atom method used in the building
        process of the Impact template.
        """

        LIGAND_PATH = get_data_file_path('ligands/malonate.pdb')
        FORCEFIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

        molecule = Molecule(LIGAND_PATH)

        openff = OpenForceField(FORCEFIELD_NAME)

        parameters = openff.parameterize(molecule, charge_method='dummy')

        topology = Topology(molecule, parameters)

        impact = Impact(topology)

        absolute_parent = impact._get_absolute_parent_atom()

        assert absolute_parent.PDB_name == '_C2_', \
            'Unexpected absolute parent atom: {}'.format(absolute_parent)

    def test_get_all_childs_of_atom(self):
        """
        It tests the _get_all_childs_of_atom method used in the building
        process of the Impact template.
        """

        LIGAND_PATH = get_data_file_path('ligands/malonate.pdb')
        FORCEFIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

        molecule = Molecule(LIGAND_PATH)

        openff = OpenForceField(FORCEFIELD_NAME)

        parameters = openff.parameterize(molecule, charge_method='dummy')

        topology = Topology(molecule, parameters)

        impact = Impact(topology)

        absolute_parent = impact._get_absolute_parent_atom()

        childs = impact._get_all_childs_of_atom(absolute_parent, 'side chain')

        assert [a.PDB_name for a in childs] == \
            ['_C1_', '_C3_'], \
            'Unexpected side-chain-child atoms: {}'.format(childs)

        childs = impact._get_all_childs_of_atom(absolute_parent, 'core')

        assert [a.PDB_name for a in childs] == \
            ['_H1_', '_H2_'], \
            'Unexpected core-child atoms: {}'.format(childs)

    def test_get_core_atoms(self):
        """
        It tests the _get_core_atoms method used in the building
        process of the Impact template.
        """

        LIGAND_PATH = get_data_file_path('ligands/malonate.pdb')
        FORCEFIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

        molecule = Molecule(LIGAND_PATH)

        openff = OpenForceField(FORCEFIELD_NAME)

        parameters = openff.parameterize(molecule, charge_method='dummy')

        topology = Topology(molecule, parameters)

        impact = Impact(topology)

        core_atoms = impact._get_core_atoms()

        assert [a.PDB_name for a in core_atoms] == \
            ['_C2_', '_H1_', '_H2_'], \
            'Unexpected core atoms: {}'.format(core_atoms)
