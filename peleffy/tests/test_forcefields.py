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
        calculator = openff._get_charge_calculator(None)(dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.Am1bccCalculator), \
            'Invalid default charge calculator: ' \
            + '{}'.format(type(calculator))

        # Check custom selection 1
        openff = OpenForceField(self.FORCE_FIELD_NAME)
        calculator = openff._get_charge_calculator('gasteiger')(dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.GasteigerCalculator), \
            'Invalid custom selection 1 for the charge calculator'

        # Check custom selection 1
        openff = OpenForceField(self.FORCE_FIELD_NAME)
        calculator = openff._get_charge_calculator('opls2005')(dummy_mol)

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
        calculator = oplsff._get_charge_calculator(None)(dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.OPLSChargeCalculator), \
            'Invalid default charge calculator: ' \
            + '{}'.format(type(calculator))

        # Check custom selection 1
        oplsff = OPLS2005ForceField()
        calculator = oplsff._get_charge_calculator('gasteiger')(dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.GasteigerCalculator), \
            'Invalid custom selection 1 for the charge calculator'

        # Check custom selection 1
        oplsff = OPLS2005ForceField()
        calculator = oplsff._get_charge_calculator('am1bcc')(dummy_mol)

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

    def test_add_solvent_parameters(self):
        """
        It tests the function that adds the solvent parameters to
        the OPLSParameters collection.
        """

        from simtk import unit
        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper

        # Using a standard atom type
        params1 = OPLS2005ParameterWrapper(
            {'atom_names': [' C1 ', ' H1 ', ' H2 ', ' H3 ', ' H4 '],
             'atom_types': ['CT', 'HC', 'HC', 'HC', 'HC'],
             'charges': [-0.24, 0.06, 0.06, 0.06, 0.06],
             'sigmas': [3.5, 2.5, 2.5, 2.5, 2.5],
             'epsilons': [0.066, 0.03, 0.03, 0.03, 0.03]})

        # Using a similar atom type
        params2 = OPLS2005ParameterWrapper(
            {'atom_names': [' C1 ', ' H1 ', ' H2 ', ' H3 ', ' H4 '],
             'atom_types': ['C3M', 'HC', 'HC', 'HC', 'HC'],
             'charges': [-0.24, 0.06, 0.06, 0.06, 0.06],
             'sigmas': [3.5, 2.5, 2.5, 2.5, 2.5],
             'epsilons': [0.066, 0.03, 0.03, 0.03, 0.03]})

        # Using a default atom type
        params3 = OPLS2005ParameterWrapper(
            {'atom_names': [' C1 ', ' H1 ', ' H2 ', ' H3 ', ' H4 '],
             'atom_types': ['XX', 'HC', 'HC', 'HC', 'HC'],
             'charges': [-0.24, 0.06, 0.06, 0.06, 0.06],
             'sigmas': [3.5, 2.5, 2.5, 2.5, 2.5],
             'epsilons': [0.066, 0.03, 0.03, 0.03, 0.03]})

        OPLS2005ParameterWrapper._add_solvent_parameters(params1)
        OPLS2005ParameterWrapper._add_solvent_parameters(params2)
        OPLS2005ParameterWrapper._add_solvent_parameters(params3)

        assert params1['SGB_radii'][0] == \
            unit.Quantity(1.975, unit.angstrom), 'Unexpected SGB radius'
        assert params1['vdW_radii'][0] == \
            unit.Quantity(1.750, unit.angstrom), 'Unexpected vdW radius'
        assert params1['gammas'][0] == 0.005000000, 'Unexpected gamma'
        assert params1['alphas'][0] == -0.741685710, 'Unexpected alpha'

        assert params2['SGB_radii'][0] == \
            unit.Quantity(2.002, unit.angstrom), 'Unexpected SGB radius'
        assert params2['vdW_radii'][0] == \
            unit.Quantity(1.775, unit.angstrom), 'Unexpected vdW radius'
        assert params2['gammas'][0] == 0.023028004, 'Unexpected gamma'
        assert params2['alphas'][0] == -0.852763146, 'Unexpected alpha'

        assert params3['SGB_radii'][0] == \
            unit.Quantity(1.500, unit.angstrom), 'Unexpected SGB radius'
        assert params3['vdW_radii'][0] == \
            unit.Quantity(1.250, unit.angstrom), 'Unexpected vdW radius'
        assert params3['gammas'][0] == 0.005000000, 'Unexpected gamma'
        assert params3['alphas'][0] == 0.000000000, 'Unexpected alpha'


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
        calculator = hybridff._get_charge_calculator(None)(dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.Am1bccCalculator), \
            'Invalid default charge calculator: ' \
            + '{}'.format(type(calculator))

        # Check custom selection 1
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        calculator = hybridff._get_charge_calculator('gasteiger')(dummy_mol)

        assert isinstance(
            calculator,
            peleffy.forcefield.calculators.GasteigerCalculator), \
            'Invalid custom selection 1 for the charge calculator'

        # Check custom selection 1
        hybridff = OpenFFOPLS2005ForceField(self.FORCE_FIELD_NAME)
        calculator = hybridff._get_charge_calculator('opls2005')(dummy_mol)

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
        hybridff.set_torsion_parameters('OPLS2005')

        check(hybridff, molecule, ffld_file,
              get_data_file_path(
                  'tests/ETL_openff-1.2.1_opls2005_parameters2.json'))
