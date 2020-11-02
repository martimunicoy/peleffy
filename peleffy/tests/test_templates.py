"""
This module contains the tests to check the peleffy's Impact template builder.
"""

import pytest

from peleffy.utils import get_data_file_path
from peleffy.topology import Molecule
from peleffy.template import Impact
from peleffy.tests.utils import compare_files


class TestImpactTemplate(object):
    """It tests the Impact class that writes a PELE's Impact template"""

    def test_input(self):
        """
        It tests that the molecule given to Impact() is of the correct
        format, peleffy.topology.Molecule and it has been parameterized
        with a forcefield.
        """
        LIGAND_PATH = 'ligands/BNZ.pdb'

        ligand_path = get_data_file_path(LIGAND_PATH)
        molecule = Molecule(ligand_path)

        # Impact() gets nothing as argument
        with pytest.raises(TypeError):
            _ = Impact()

        # Impact() gets a non topology.Molecule object
        with pytest.raises(TypeError):
            _ = Impact('passing a str instead of a Molecule')

        # The molecule is not parameterized
        with pytest.raises(Exception):
            _ = Impact(molecule)

    def _prepare_molecule_OPLS(self, pdb_name, ffld_name,
                               molecule_tag=None):
        """
        It initialices and parametrizes a molecule using OPLS.

        Parameters
        ----------
        pdb_name: str
            The relative path to the PDB of the molecule.
        ffld_name : str
            The relative path to the TXT file containing OPLS parameters.
        molecule_tag : str
            The tag to set to the molecule. Default is None

        """
        from peleffy.forcefield import (OPLS2005ForceField,
                                        OPLS2005ParameterWrapper)

        FORCE_FIELD_NAME = 'OPLS2005'
        # Loads the molecule and the force field
        pdb_path = get_data_file_path(pdb_name)
        if molecule_tag is not None:
            molecule = Molecule(pdb_path, tag=molecule_tag)
        else:
            molecule = Molecule(pdb_path)
        oplsff = OPLS2005ForceField(FORCE_FIELD_NAME)

        # Workaround to skip Schrodinger dependency
        ffld_file = get_data_file_path(ffld_name)
        with open(ffld_file) as f:
            ffld_output = f.read()
        oplsff._parameters = \
            OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                      ffld_output)

        # Set force field and parameterize molecule
        molecule.set_forcefield(oplsff)
        molecule.parameterize()

        return molecule

    def test_writer_OFF(self):
        """
        It tests the writer attribute of the Impact class using OFF
        to parameterize.
        """

        TEMPLATE_METZ = get_data_file_path('tests/metz')
        TEMPLATE_MATZ = get_data_file_path('tests/malz')
        TEMPLATE_ETLZ = get_data_file_path('tests/etlz')

        # Generates the template for methane
        pdb_path = get_data_file_path('ligands/methane.pdb')
        m = Molecule(pdb_path)
        m.parameterize('openff_unconstrained-1.2.1.offxml')
        impact = Impact(m)
        impact.write('metz')

        # Compare the reference template and the generated template
        compare_files(file1=TEMPLATE_METZ, file2='metz')

        # Generates the template for malonate
        pdb_path = get_data_file_path('ligands/malonate.pdb')
        m = Molecule(pdb_path)
        m.parameterize('openff_unconstrained-1.2.1.offxml')
        impact = Impact(m)
        impact.write('malz')

        # Compare the reference template and the generated template
        compare_files(file1=TEMPLATE_MATZ, file2='malz')

        # Generates the template for ethylene
        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        m = Molecule(pdb_path, tag='ETL')  # Note that in this case we are assigning a tag to the molecule which will be used in the Impact template
        m.parameterize('openff_unconstrained-1.2.1.offxml')
        impact = Impact(m)
        impact.write('etlz')

        # Compare the reference template and the generated template
        compare_files(file1=TEMPLATE_ETLZ, file2='etlz')

    def test_writer_OPLS(self):
        """
        It tests the writer attribute of the Impact class using OPLS to parameterize.
        """
        TEMPLATE_METZ_OPLS = get_data_file_path('tests/OPLS_metz')
        TEMPLATE_MALZ_OPLS = get_data_file_path('tests/OPLS_malz')
        TEMPLATE_ETLZ_OPLS = get_data_file_path('tests/OPLS_etlz')

        # Generates the template for methane using OPLS
        m = self._prepare_molecule_OPLS(pdb_name='ligands/methane.pdb',
                                        ffld_name='tests/MET_ffld_output.txt')
        impact = Impact(m)
        impact.write('metz')

        # Compare the reference template and the generated template
        compare_files(file1=TEMPLATE_METZ_OPLS, file2='metz')

        # Generates the template for malonate
        m = self._prepare_molecule_OPLS(pdb_name='ligands/malonate.pdb',
                                        ffld_name='tests/MAL_ffld_output.txt')
        impact = Impact(m)
        impact.write('malz')

        # Compare the reference template and the generated template
        compare_files(file1=TEMPLATE_MALZ_OPLS, file2='malz')

        # Generates the template for ethylene
        m = self._prepare_molecule_OPLS(pdb_name='ligands/ethylene.pdb',
                                        ffld_name='tests/ETL_ffld_output.txt',
                                        molecule_tag='ETL')
        impact = Impact(m)
        impact.write('etlz')

        # Compare the reference template and the generated template
        compare_files(file1=TEMPLATE_ETLZ_OPLS, file2='etlz')
