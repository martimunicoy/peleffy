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
		TEMPLATE_METZ = '/home/lauramalo/repos/offpele/offpele/tests/reference_templates/metz'
		TEMPLATE_MATZ = '/home/lauramalo/repos/offpele/offpele/tests/reference_templates/malz'
		TEMPLATE_ETLZ = '/home/lauramalo/repos/offpele/offpele/tests/reference_templates/etlz'

		# Generates the tempalte for methane
		pdb_path = get_data_file_path('ligands/methane.pdb')
		m = Molecule(pdb_path)
		m.parameterize('openff_unconstrained-1.2.1.offxml')
		impact = Impact(m)
		impact.write('metz')

		# Compare the reference template and the generated template
		with open(TEMPLATE_METZ, 'r') as f1: 
			dataA = f1.readlines()[3:]
		with open('metz', 'r') as f2: 
			dataB = f2.readlines()[3:]
		assert dataA == dataB

		# Generates the template for malonate
		pdb_path = get_data_file_path('ligands/malonate.pdb')
		m = Molecule(pdb_path)
		m.parameterize('openff_unconstrained-1.2.1.offxml')
		impact = Impact(m)
		impact.write('malz')

		# Compare the reference template and the generated template
		with open(TEMPLATE_MATZ, 'r') as f1: 
			dataA = f1.readlines()[3:]
		with open('malz', 'r') as f2: 
			dataB = f2.readlines()[3:]
		assert dataA == dataB

		# Generates the template for ethylene
		pdb_path = get_data_file_path('ligands/ethylene.pdb')
		m = Molecule(pdb_path, tag='ETL')  # Note that in this case we are assigning a tag to the molecule which will be used in the Impact template
		m.parameterize('openff_unconstrained-1.2.1.offxml')
		impact = Impact(m)
		impact.write('etlz')

		# Compare the reference template and the generated template
		with open(TEMPLATE_ETLZ, 'r') as f1: 
			dataA = f1.readlines()[3:]
		with open('etlz', 'r') as f2: 
			dataB = f2.readlines()[3:]
		assert dataA == dataB


	










