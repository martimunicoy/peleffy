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


	def _compare_files(self, FILE1, FILE2): 
		"""
			It comapares two files line by line and gives an AssertionError if there 
			is any difference. 
		"""
		with open(FILE1, 'r') as f1: 
			dataA = f1.readlines()[3:]
		with open(FILE2, 'r') as f2: 
			dataB = f2.readlines()[3:]
		linesA = str(dataA).split('\n')
		linesB = str(dataB).split('\n')
		assert len(linesA) == len(linesB), 'Number of lines do not match'
		for i, (lineA, lineB) in enumerate(zip(linesA, linesB)):
			assert lineA == lineB, 'Found different lines at line {}:' + lineA + ' - ' + lineB

	def test_writer(self):
		"""
		It tests the writer attribute of the Impact class. 
		"""

		FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'

		TEMPLATE_METZ = get_data_file_path('tests/metz')
		TEMPLATE_MATZ = get_data_file_path('tests/malz')
		TEMPLATE_ETLZ = get_data_file_path('tests/etlz')

		# Generates the tempalte for methane
		pdb_path = get_data_file_path('ligands/methane.pdb')
		m = Molecule(pdb_path)
		m.parameterize('openff_unconstrained-1.2.1.offxml')
		impact = Impact(m)
		impact.write('metz')

		# Compare the reference template and the generated template
		self._compare_files(FILE1 = TEMPLATE_METZ, FILE2 = 'metz')

		# Generates the template for malonate
		pdb_path = get_data_file_path('ligands/malonate.pdb')
		m = Molecule(pdb_path)
		m.parameterize('openff_unconstrained-1.2.1.offxml')
		impact = Impact(m)
		impact.write('malz')

		# Compare the reference template and the generated template
		self._compare_files(FILE1 = TEMPLATE_MATZ, FILE2 = 'malz')

		# Generates the template for ethylene
		pdb_path = get_data_file_path('ligands/ethylene.pdb')
		m = Molecule(pdb_path, tag='ETL')  # Note that in this case we are assigning a tag to the molecule which will be used in the Impact template
		m.parameterize('openff_unconstrained-1.2.1.offxml')
		impact = Impact(m)
		impact.write('etlz')

		# Compare the reference template and the generated template
		self._compare_files(FILE1 = TEMPLATE_ETLZ, FILE2 = 'etlz')







