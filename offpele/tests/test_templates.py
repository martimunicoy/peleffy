"""
This module contains the tests to check the offpele's Impact template builder.
"""

import pytest


from offpele.utils import get_data_file_path
from offpele.topology import Molecule
from offpele.template import Impact

class TestImpactTemplate(object): 

	"""It tests the Impact class that writes a PELE's Impact template"""

	def test_input(self): 
		"""
		It tests that the molecule given to Impact() is of the correct format, offpele.topology.Molecule
		and it has been parameterized with a forcefield. 
		"""
		
		LIGAND_PATH = 'ligands/BNZ.pdb'
		FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'

		ligand_path = get_data_file_path(LIGAND_PATH)
		molecule = Molecule(ligand_path)
		
		# The molecule is not parameterized
		with pytest.raises(AssertionError): 
			impact = Impact(molecule)

		# Impact() gets nothing as argument
		molecule.parameterize(FORCEFILED_NAME)
		with pytest.raises(Exception): 
			impact = Impact()

		# Molecule is not a topology.molecule object
		molecule_string= Molecule('non_molecule_type')
		molecule_string.parameterize(FORCEFIELD_NAME)
		with pytest.raises(TypeError):
			impact = Impact(molecule_string)


	def test_writer(self):
		"""
		It tests the writer attribute of the Impact class. 
		"""

		LIGAND_PATH = 'ligands/BNZ.pdb'
		FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'
		TEMPLATE_EXAMPLE = 'ghfdh'

		ligand_path = get_data_file_path(LIGAND_PATH)
		molecule = Molecule(ligand_path)
		molecule.parameterize(FORCEFIELD_NAME)
		impact = Impact(molecule)
		impact.write('molz')

		with open(TEMPLATE_EXAMPLE, 'r') as f1: 
			contentA = set(f1)
		with open('molz', 'r') as f2: 
			contentB = set(f2)
		assert contentA == contentB












