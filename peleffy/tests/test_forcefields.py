"""
This module contains the tests to check all available force fields in
peleffy.
"""


class TestOpenForceField(object):
    """
    It wraps all tests that check the OpenForceField class.
    """

    FORCE_FIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

    def test_name(self):
        """It checks the name assignment."""

        from peleffy.forcefield import OpenForceField

        openff = OpenForceField(self.FORCE_FIELD_NAME)

        assert openff.name == self.FORCE_FIELD_NAME, \
            'Unexpected force field name'

    def test_type(self):
        """It checks the type assignment."""

        from peleffy.forcefield import OpenForceField

        openff = OpenForceField(self.FORCE_FIELD_NAME)

        assert openff.type == 'OpenFF', \
            'Unexpected force field type'

    def test_charge_calculator_selector(self):
        """It checks the charge calculator selector."""
        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenForceField
        import peleffy

        dummy_mol = Molecule()

        # Check default selection
        openff = OpenForceField(self.FORCE_FIELD_NAME)
        calculator = openff._get_charge_calculator(None, dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.Am1bccCalculator), \
            'Invalid default charge calculator: ' \
            + '{}'.format(type(calculator))

        # Check custom selection 1
        openff = OpenForceField(self.FORCE_FIELD_NAME)
        calculator = openff._get_charge_calculator('gasteiger', dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.GasteigerCalculator), \
            'Invalid custom selection 1 for the charge calculator'

        # Check custom selection 1
        openff = OpenForceField(self.FORCE_FIELD_NAME)
        calculator = openff._get_charge_calculator('opls2005', dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.OPLSChargeCalculator), \
            'Invalid custom selection 2 for the charge calculator'

    def test_parameterizer(self):
        """It checks the parameterized method."""

        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenForceField
        from peleffy.utils import (get_data_file_path,
                                   convert_all_quantities_to_string)
        from .utils import compare_dicts
        import json

        # Load molecule 1
        molecule = Molecule(get_data_file_path('ligands/methane.pdb'))
        openff = OpenForceField(self.FORCE_FIELD_NAME)

        # Obtain force field parameters
        parameters = openff.parameterize(molecule)

        writable_parameters = convert_all_quantities_to_string(parameters)

        reference_file = get_data_file_path(
            'tests/MET_openff-1.2.1_parameters.json')

        with open(reference_file) as f:
            compare_dicts(writable_parameters, json.load(f))

        # Load molecule
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))
        openff = OpenForceField(self.FORCE_FIELD_NAME)

        # Obtain force field parameters
        parameters = openff.parameterize(molecule)

        writable_parameters = convert_all_quantities_to_string(parameters)

        reference_file = get_data_file_path(
            'tests/ETL_openff-1.2.1_parameters.json')

        with open(reference_file) as f:
            compare_dicts(writable_parameters, json.load(f))

    def test_OPLS_charge_method_with_openff(self):
        """
        It tests the obtaining of OPLS2005 charges using the OpenFF
        force field.
        """
        import peleffy
        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenForceField
        from peleffy.utils import get_data_file_path

        # Load molecule
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))

        # Load OpenFF force field with the OPLS2005 charge method
        openff = OpenForceField(self.FORCE_FIELD_NAME)

        charge_calculator = openff._get_charge_calculator(
            charge_method='OPLS2005',
            molecule=molecule)

        assert isinstance(
            charge_calculator,
            peleffy.forcefield.calculators.OPLSChargeCalculator), \
            'Unexpected charge calculator method'


class TestOPLS2005ForceField(object):
    """
    It wraps all tests that check the OPLS2005ForceField class.
    """

    FORCE_FIELD_NAME = 'OPLS2005'

    def test_name(self):
        """It checks the name assignment."""

        from peleffy.forcefield import OPLS2005ForceField

        oplsff = OPLS2005ForceField()

        assert oplsff.name == self.FORCE_FIELD_NAME, \
            'Unexpected force field name'

    def test_type(self):
        """It checks the type assignment."""

        from peleffy.forcefield import OPLS2005ForceField

        oplsff = OPLS2005ForceField()

        assert oplsff.type == 'OPLS2005', \
            'Unexpected force field type'

    def test_charge_calculator_selector(self):
        """It checks the charge calculator selector."""
        from peleffy.topology import Molecule
        from peleffy.forcefield import OPLS2005ForceField
        import peleffy

        dummy_mol = Molecule()

        # Check default selection
        oplsff = OPLS2005ForceField()
        calculator = oplsff._get_charge_calculator(None, dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.OPLSChargeCalculator), \
            'Invalid default charge calculator: ' \
            + '{}'.format(type(calculator))

        # Check custom selection 1
        oplsff = OPLS2005ForceField()
        calculator = oplsff._get_charge_calculator('gasteiger', dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.GasteigerCalculator), \
            'Invalid custom selection 1 for the charge calculator'

        # Check custom selection 1
        oplsff = OPLS2005ForceField()
        calculator = oplsff._get_charge_calculator('am1bcc', dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.Am1bccCalculator), \
            'Invalid custom selection 2 for the charge calculator'

    def test_parameterizer(self):
        """It checks the parameterized method."""

        from peleffy.topology import Molecule
        from peleffy.forcefield import OPLS2005ForceField
        from peleffy.utils import (get_data_file_path,
                                   convert_all_quantities_to_string)
        from .utils import compare_dicts, parameterize_opls2005
        import json

        # Load molecule 1
        molecule = Molecule(get_data_file_path('ligands/methane.pdb'))
        oplsff = OPLS2005ForceField()
        ffld_file = get_data_file_path('tests/MET_ffld_output.txt')

        parameters = parameterize_opls2005(oplsff, molecule, ffld_file)

        writable_parameters = convert_all_quantities_to_string(parameters)

        reference_file = get_data_file_path(
            'tests/MET_opls2005_parameters.json')

        with open(reference_file) as f:
            compare_dicts(writable_parameters, json.load(f))

        # Load molecule 1
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))
        oplsff = OPLS2005ForceField()
        ffld_file = get_data_file_path('tests/ETL_ffld_output.txt')

        parameters = parameterize_opls2005(oplsff, molecule, ffld_file)

        writable_parameters = convert_all_quantities_to_string(parameters)

        reference_file = get_data_file_path(
            'tests/ETL_opls2005_parameters.json')

        with open(reference_file) as f:
            compare_dicts(writable_parameters, json.load(f))


class TestOpenFFOPLS2005ForceField(object):
    """
    It wraps all tests that check the OpenFFOPLS2005ForceField class.
    """

    FORCE_FIELD_NAME = 'openff_unconstrained-1.2.1.offxml'

    def test_name(self):
        """It checks the name assignment."""

        from peleffy.forcefield import OpenFFOPLS2005ForceField

        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)

        assert hybridff.name == self.FORCE_FIELD_NAME + ' + OPLS2005', \
            'Unexpected force field name'

    def test_type(self):
        """It checks the type assignment."""

        from peleffy.forcefield import OpenFFOPLS2005ForceField

        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)

        assert hybridff.type == 'OpenFF + OPLS2005', \
            'Unexpected force field type'

    def test_charge_calculator_selector(self):
        """It checks the charge calculator selector."""
        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenFFOPLS2005ForceField
        import peleffy

        dummy_mol = Molecule()

        # Check default selection
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        calculator = hybridff._get_charge_calculator(None, dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.Am1bccCalculator), \
            'Invalid default charge calculator: ' \
            + '{}'.format(type(calculator))

        # Check custom selection 1
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        calculator = hybridff._get_charge_calculator('gasteiger', dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.GasteigerCalculator), \
            'Invalid custom selection 1 for the charge calculator'

        # Check custom selection 1
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        calculator = hybridff._get_charge_calculator('opls2005', dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.OPLSChargeCalculator), \
            'Invalid custom selection 2 for the charge calculator'

    def test_parameterizer(self):
        """It checks the parameterized method."""

        def check(hybridff, molecule, ffld_file, reference_file):
            """
            It checks the parameters obtained by the force field.

            Parameters
            ----------
            hybridff : an OpenFFOPLS2005ForceField object
                The hybrid force field to employ in the parameterization
                along with the ffld_file
            molecule : a peleffy.topology.Molecule
                The peleffy's Molecule object to parameterize with the
                ffld file
            ffld_file : str
                The path to the precomputed ffld file from where the
                parameters will be extracted
            reference_file : str
                The path to the file containing the reference parameters
            """
            parameters = parameterize_openffopls2005(hybridff,
                                                     molecule,
                                                     ffld_file)

            writable_parameters = convert_all_quantities_to_string(parameters)

            with open(reference_file) as f:
                compare_dicts(writable_parameters, json.load(f))

        from peleffy.topology import Molecule
        from peleffy.forcefield import OpenFFOPLS2005ForceField
        from peleffy.utils import (get_data_file_path,
                                   convert_all_quantities_to_string)
        from .utils import compare_dicts, parameterize_openffopls2005
        import json

        # Load molecule 1
        molecule = Molecule(get_data_file_path('ligands/methane.pdb'))
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        ffld_file = get_data_file_path('tests/MET_ffld_output.txt')

        # 1st check
        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/MET_openff-1.2.1_opls2005_parameters1.json'))

        # 2nd check
        hybridff.set_nonbonding_parameters('OPLS2005')

        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/MET_openff-1.2.1_opls2005_parameters2.json'))

        # 3rd check
        hybridff.set_nonbonding_parameters('OpenFF')
        hybridff.set_bond_parameters('OPLS2005')

        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/MET_openff-1.2.1_opls2005_parameters3.json'))

        # 4th check
        hybridff.set_nonbonding_parameters('OpenFF')
        hybridff.set_bond_parameters('OpenFF')
        hybridff.set_angle_parameters('OPLS2005')

        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/MET_openff-1.2.1_opls2005_parameters4.json'))

        # 5th check
        hybridff.set_nonbonding_parameters('OPLS2005')
        hybridff.set_bond_parameters('OPLS2005')
        hybridff.set_angle_parameters('OPLS2005')

        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/MET_openff-1.2.1_opls2005_parameters5.json'))

        # Load molecule 2
        molecule = Molecule(get_data_file_path('ligands/ethylene.pdb'))
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        ffld_file = get_data_file_path('tests/ETL_ffld_output.txt')

        # 1st check
        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/ETL_openff-1.2.1_opls2005_parameters1.json'))

        # 2nd check
        hybridff.set_nonbonding_parameters('OpenFF')
        hybridff.set_dihedral_parameters('OPLS2005')

        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/ETL_openff-1.2.1_opls2005_parameters2.json'))


class TestForceFieldSelector(object):
    """
    It wraps all tests that check the force field selector class.
    """

    def test_get_by_name(self):
        """It checks the get_by_name method."""

        from peleffy.forcefield import ForceFieldSelector
        from peleffy.forcefield import OpenForceField
        from peleffy.forcefield import OPLS2005ForceField

        selector = ForceFieldSelector()

        forcefield = selector.get_by_name('openff_unconstrained-1.0.0.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('openff_unconstrained-1.0.1.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('openff_unconstrained-1.1.0.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('openff_unconstrained-1.1.1.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('openff_unconstrained-1.2.0.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('openff_unconstrained-1.2.1.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('openff_unconstrained-1.3.0.offxml')

        assert isinstance(forcefield, OpenForceField), \
            'Unexpected force field type'

        forcefield = selector.get_by_name('OPLS2005')

        assert isinstance(forcefield, OPLS2005ForceField), \
            'Unexpected force field type'
