"""
This module contains the tests to check peleffy's parameters.
"""

import pytest

import tempfile
from simtk import unit
import numpy as np

from .utils import (SET_OF_LIGAND_PATHS, apply_PELE_dihedral_equation,
                    apply_OFF_dihedral_equation, check_CHO_charges,
                    compare_files)
from peleffy.topology import Molecule, Bond, Angle, Proper
from peleffy.template.impact import (WritableBond, WritableAngle,
                                     WritableProper, WritableImproper)
from peleffy.forcefield import OpenForceField, OPLS2005ForceField
from peleffy.forcefield.parameters import BaseParameterWrapper
from peleffy.topology import Topology
from peleffy.template import Impact
from peleffy.utils import get_data_file_path, temporary_cd


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'


class TestWrapper(object):
    """
    It contains all the tests that validate the parameter wrapper class.
    """

    def test_forcefield_name_assignment(self):
        """
        It validates the force field name assignment of the parameter
        wrapper.
        """
        from peleffy.forcefield.parameters \
            import (BaseParameterWrapper, OpenForceFieldParameterWrapper,
                    OPLS2005ParameterWrapper, OpenFFOPLS2005ParameterWrapper)

        p = BaseParameterWrapper()

        assert p.forcefield_name == '', \
            'Unexpected force field name found in the parameters wrapper'

        p = BaseParameterWrapper(
            forcefield_name='openff_unconstrained-1.2.1.offxml')

        assert p.forcefield_name == 'openff_unconstrained-1.2.1.offxml', \
            'Unexpected force field name found in the parameters wrapper'

        p = OpenForceFieldParameterWrapper()

        assert p.forcefield_name == 'OpenFF', \
            'Unexpected force field name found in the parameters wrapper'

        p = OpenForceFieldParameterWrapper(
            forcefield_name='openff_unconstrained-1.2.1.offxml')

        assert p.forcefield_name == 'openff_unconstrained-1.2.1.offxml', \
            'Unexpected force field name found in the parameters wrapper'

        p = OPLS2005ParameterWrapper()

        assert p.forcefield_name == 'OPLS2005', \
            'Unexpected force field name found in the parameters wrapper'

        p = OpenForceFieldParameterWrapper(forcefield_name='opls')

        assert p.forcefield_name == 'opls', \
            'Unexpected force field name found in the parameters wrapper'

        p = OpenFFOPLS2005ParameterWrapper()

        assert p.forcefield_name == 'Openff + OPLS2005', \
            'Unexpected force field name found in the parameters wrapper'

        p = OpenForceFieldParameterWrapper(
            forcefield_name='openff_unconstrained-1.2.1.offxml')

        assert p.forcefield_name == 'openff_unconstrained-1.2.1.offxml', \
            'Unexpected force field name found in the parameters wrapper'

    def test_comparison(self):
        """It tests the comparison between parameter wrapper."""
        from peleffy.forcefield.parameters import BaseParameterWrapper

        p = BaseParameterWrapper()
        p2 = BaseParameterWrapper()
        p3 = BaseParameterWrapper(
            forcefield_name='openff_unconstrained-1.2.1.offxml')
        p4 = BaseParameterWrapper({'atom_names': ' C1 '})
        p5 = BaseParameterWrapper({'atom_names': ' C1 '})
        p6 = BaseParameterWrapper(
            {'atom_names': ' C1 '},
            forcefield_name='openff_unconstrained-1.2.1.offxml')

        assert p == p2, 'Unexpected inequality between \'p\' and \'p2\''
        assert p != p3, 'Unexpected equality between \'p\' and \'p3\''
        assert p != p4, 'Unexpected equality between \'p\' and \'p4\''
        assert p4 == p5, 'Unexpected inequality between \'p4\' and \'p5\''
        assert p4 != p6, 'Unexpected equality between \'p4\' and \'p6\''

    def test_to_string(self):
        """It tests the string representation of the parameter wrapper."""

        pdb_path = get_data_file_path('ligands/methane.pdb')
        molecule = Molecule(pdb_path)
        ff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        parameters = ff.parameterize(molecule,
                                     charge_method='gasteiger')
        p_string = parameters.to_string()

        ref_string_file = get_data_file_path(
            'tests/MET_parameters_to_string.txt')
        with open(ref_string_file) as f:
            ref_string = f.read().strip('\n')

        p_lines = p_string.split('\n')
        ref_lines = ref_string.split('\n')

        assert len(p_lines) == len(ref_lines), \
            'Unexpected number of lines: ' \
            + str(len(p_lines)) + ', expected ' + str(len(ref_lines))

        for p_line, ref_line in zip(p_lines, ref_lines):
            assert p_line == ref_line, \
                'Unexpected string representation found in line: ' \
                + '\'{}\' '.format(p_line) \
                + 'which does not match with \'{}\''.format(ref_line)

    def test_to_json(self):
        """It tests the json representation of the parameter wrapper."""

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                pdb_path = get_data_file_path('ligands/methane.pdb')
                molecule = Molecule(pdb_path)
                ff = OpenForceField('openff_unconstrained-1.2.1.offxml')
                parameters = ff.parameterize(molecule,
                                             charge_method='gasteiger')
                parameters.to_json('parameters_to_check.json')

                ref_json_file = get_data_file_path(
                    'tests/MET_parameters_to_json.json')
                compare_files('parameters_to_check.json', ref_json_file)

    def test_from_json(self):
        """It tests the loading function from json files into the parameter
        wrapper"""

        def compare_BaseParameterWrapper(params1, params2):
            """
            It compares two BaseParameterWrapper objects.
            """
            import numpy as np
            for p in params1.keys():
                if p == 'charges':      #charges have to be handle differently
                    assert params1[p].__eq__(params2[p]).all() == \
                            np.full((len(params1[p])), True, dtype=bool).all()
                else:
                    assert params1[p] == params2[p]

        with tempfile.TemporaryDirectory() as tmpdir:
            with temporary_cd(tmpdir):
                #Force Field to parameterize the molecules
                ff = OpenForceField('openff_unconstrained-1.2.1.offxml')

                # Test for malonate
                m = Molecule(get_data_file_path('ligands/malonate.pdb'))
                params_ref = ff.parameterize(m)
                params_ref.to_json('parameters_mal.json')

                from peleffy.forcefield.parameters import \
                                OpenForceFieldParameterWrapper
                wrapper_off = OpenForceFieldParameterWrapper()
                params_load = wrapper_off.from_json('parameters_mal.json')

                compare_BaseParameterWrapper(params_ref, params_load)

                # Test for methane
                m = Molecule(get_data_file_path('ligands/methane.pdb'))
                params_ref = ff.parameterize(m)
                params_ref.to_json('parameters_met.json')

                from peleffy.forcefield.parameters import \
                                OpenForceFieldParameterWrapper
                wrapper_off = OpenForceFieldParameterWrapper()
                params_load = wrapper_off.from_json('parameters_met.json')

                compare_BaseParameterWrapper(params_ref, params_load)

                # Test for ethylene
                m = Molecule(get_data_file_path('ligands/ethylene.pdb'))
                params_ref = ff.parameterize(m)
                params_ref.to_json('parameters_etl.json')

                from peleffy.forcefield.parameters import \
                                OpenForceFieldParameterWrapper
                wrapper_off = OpenForceFieldParameterWrapper()
                params_load = wrapper_off.from_json('parameters_etl.json')

                compare_BaseParameterWrapper(params_ref, params_load)

    def test_from_impact_template(self):
        """
        It tests the method to generate a parameter wrapper out of an
        impact template.
        """
        def test_generate_OpenForceFieldParameterWrapper(molecule,
                                                         impact_template_path):
            """
            It tests the method to return a OpenForceFieldParameterWrapper
            object from an Impact file by comparing the reference Impact file
            and the generated from the obtained Wrapper.
            """

            from peleffy.forcefield.parameters import \
                OpenForceFieldParameterWrapper

            with tempfile.TemporaryDirectory() as tmpdir:
                with temporary_cd(tmpdir):

                    # Assign parameters to BaseParametersWrapper object from
                    # an Impact Template
                    wrapper_off = OpenForceFieldParameterWrapper()
                    parameters = wrapper_off.from_impact_template(
                        molecule, impact_template_path)

                    # Generate the new impact template
                    topology = Topology(molecule, parameters)
                    impact = Impact(topology)
                    impact.to_file('unlz_generated')

                    # Compare the reference and generated Impact templates
                    compare_files(impact_template_path, 'unlz_generated')

        def test_generate_OPLS2005ParameterWrapper(molecule,
                                                   impact_template_path):
            """
            It tests the method to return a OPLS2005FieldParameterWrapper
            object from an Impact file by comparing the reference Impact file
            and the generated from the obtained Wrapper.
            """

            from peleffy.forcefield.parameters import OPLS2005ParameterWrapper

            with tempfile.TemporaryDirectory() as tmpdir:
                with temporary_cd(tmpdir):

                    # Assign parameters to BaseParametersWrapper object from
                    # an Impact Template
                    wrapper_opls = OPLS2005ParameterWrapper()
                    parameters = wrapper_opls.from_impact_template(
                        molecule, impact_template_path)

                    # Generate the new impact template
                    topology = Topology(molecule, parameters)
                    impact = Impact(topology)
                    impact.to_file('unlz_generated')

                    # Compare the reference and generated Impact templates
                    compare_files(impact_template_path, 'unlz_generated')

        # Test with OFF parametrization (for ethylene)
        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        molecule = Molecule(pdb_path, tag='ETL')
        impact_template_path = get_data_file_path('tests/etlz')
        test_generate_OpenForceFieldParameterWrapper(molecule,
                                                     impact_template_path)

        # Test with OFF parametrization (for methane)
        pdb_path = get_data_file_path('ligands/methane.pdb')
        molecule = Molecule(pdb_path)
        impact_template_path = get_data_file_path('tests/metz')
        test_generate_OpenForceFieldParameterWrapper(molecule,
                                                     impact_template_path)

        # Test with OFF parametrization (for malonate)
        pdb_path = get_data_file_path('ligands/malonate.pdb')
        molecule = Molecule(pdb_path)
        impact_template_path = get_data_file_path('tests/malz')
        test_generate_OpenForceFieldParameterWrapper(molecule,
                                                     impact_template_path)

        # Test with OPLS parametrization (for ethylene)
        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        molecule = Molecule(pdb_path, tag='ETL')
        impact_template_path = get_data_file_path('tests/OPLS_etlz')
        test_generate_OPLS2005ParameterWrapper(molecule,
                                               impact_template_path)

        # The molecule and Impact template do no represent the same chemical
        # entity
        with pytest.raises(ValueError):
            pdb_path = get_data_file_path('ligands/ethylene.pdb')
            molecule = Molecule(pdb_path, tag='ETL')
            impact_template_path = get_data_file_path('tests/metz')
            test_generate_OpenForceFieldParameterWrapper(molecule,
                                                         impact_template_path)


class TestBonds(object):
    """
    It wraps all tests that involve bond parameters.
    """

    def test_OFF_parameters(self):
        """
        It checks the standard Open Force Field parameterization for
        bonds.
        """

        MAX_THRESHOLD = 1e-3

        # Load benzene ring
        molecule = Molecule(smiles='c1ccccc1')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Check resulting parameters
        expected_lengths = {(0, 1): 1.389761064076,
                            (0, 5): 1.389761064076,
                            (0, 6): 1.085503378387,
                            (1, 2): 1.389761064076,
                            (1, 7): 1.085503378387,
                            (2, 3): 1.389761064076,
                            (2, 8): 1.085503378387,
                            (3, 4): 1.389761064076,
                            (3, 9): 1.085503378387,
                            (4, 5): 1.389761064076,
                            (4, 10): 1.085503378387,
                            (5, 11): 1.085503378387}

        expected_ks = {(0, 1): 672.0064645264,
                       (0, 5): 672.0064645264,
                       (0, 6): 808.4160937,
                       (1, 2): 672.0064645264,
                       (1, 7): 808.4160937,
                       (2, 3): 672.0064645264,
                       (2, 8): 808.4160937,
                       (3, 4): 672.0064645264,
                       (3, 9): 808.4160937,
                       (4, 5): 672.0064645264,
                       (4, 10): 808.4160937,
                       (5, 11): 808.4160937}

        for bond in parameters['bonds']:
            indexes = (bond['atom1_idx'], bond['atom2_idx'])
            expected_length = unit.Quantity(expected_lengths[indexes],
                                            unit.angstrom)
            expected_k = unit.Quantity(expected_ks[indexes],
                                       unit.kilocalorie
                                       / (unit.angstrom ** 2 * unit.mole))

            assert bond['eq_dist'] - expected_length \
                < unit.Quantity(MAX_THRESHOLD, unit.angstrom), \
                'Invalid length for bond {}'.format(bond)
            assert bond['spring_constant'] - expected_k \
                < unit.Quantity(MAX_THRESHOLD, unit.kilocalorie
                                / (unit.angstrom ** 2 * unit.mole)), \
                'Invalid k for bond {}'.format(bond)

    def test_Impact_writable_parameters(self):
        """
        It checks the Impact writable representation for bonds.
        """

        # Load benzene ring
        molecule = Molecule(smiles='CC=O')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Generate topology
        topology = Topology(molecule, parameters)

        expected_parameters = list([[1, 2, 332.5750972667, 1.523640340452],
                                    [1, 4, 376.8940758588, 1.094223427522],
                                    [1, 5, 376.8940758588, 1.094223427522],
                                    [1, 6, 376.8940758588, 1.094223427522],
                                    [2, 3, 608.3286693405, 1.225108345696],
                                    [2, 7, 404.20804685, 1.085503378387]])

        # Check resulting parameters
        for bond in topology.bonds:
            w_bond = WritableBond(bond)
            w_parameters = [attr[1] for attr in list(w_bond)]
            assert w_parameters in expected_parameters, \
                'Invalid writable bond parameters {}'.format(w_parameters)


class TestAngles(object):
    """
    It wraps all tests that involve angle parameters.
    """

    def test_OFF_parameters(self):
        """
        It checks the standard Open Force Field parameterization for
        angles.
        """

        MAX_THRESHOLD = 1e-3

        # Load benzene ring
        molecule = Molecule(smiles='c1ccccc1')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Check resulting parameters
        expected_angles = {(0, 1, 2): 128.2771922378,
                           (0, 1, 7): 133.1339832262,
                           (0, 5, 4): 128.2771922378,
                           (0, 5, 11): 133.1339832262,
                           (1, 0, 5): 128.2771922378,
                           (1, 0, 6): 133.1339832262,
                           (1, 2, 3): 128.2771922378,
                           (1, 2, 8): 133.1339832262,
                           (2, 1, 7): 133.1339832262,
                           (2, 3, 4): 128.2771922378,
                           (2, 3, 9): 133.1339832262,
                           (3, 2, 8): 133.1339832262,
                           (3, 4, 5): 128.2771922378,
                           (3, 4, 10): 133.1339832262,
                           (4, 3, 9): 133.1339832262,
                           (4, 5, 11): 133.1339832262,
                           (5, 0, 6): 133.1339832262,
                           (5, 4, 10): 133.1339832262}

        expected_ks = {(0, 1, 2): 157.3576298529,
                       (0, 1, 7): 68.40592742547,
                       (0, 5, 4): 157.3576298529,
                       (0, 5, 11): 68.40592742547,
                       (1, 0, 5): 157.3576298529,
                       (1, 0, 6): 68.40592742547,
                       (1, 2, 3): 157.3576298529,
                       (1, 2, 8): 68.40592742547,
                       (2, 1, 7): 68.40592742547,
                       (2, 3, 4): 157.3576298529,
                       (2, 3, 9): 68.40592742547,
                       (3, 2, 8): 68.40592742547,
                       (3, 4, 5): 157.3576298529,
                       (3, 4, 10): 68.40592742547,
                       (4, 3, 9): 68.40592742547,
                       (4, 5, 11): 68.40592742547,
                       (5, 0, 6): 68.40592742547,
                       (5, 4, 10): 68.40592742547}

        for angle in parameters['angles']:
            indexes = (angle['atom1_idx'], angle['atom2_idx'],
                       angle['atom3_idx'])
            expected_angle = unit.Quantity(expected_angles[indexes],
                                           unit.degree)
            expected_k = unit.Quantity(expected_ks[indexes],
                                       unit.kilocalorie
                                       / (unit.radian ** 2 * unit.mole))

            assert angle['eq_angle'] - expected_angle \
                < unit.Quantity(MAX_THRESHOLD, unit.degree), \
                'Invalid length for angle {}'.format(angle)
            assert angle['spring_constant'] - expected_k \
                < unit.Quantity(MAX_THRESHOLD, unit.kilocalorie
                                / (unit.radian ** 2 * unit.mole)), \
                'Invalid k for angle {}'.format(angle)

    def test_Impact_writable_parameters(self):
        """
        It checks the Impact writable representation for angles.
        """

        # Load benzene ring
        molecule = Molecule(smiles='CC=O')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Generate topology
        topology = Topology(molecule, parameters)

        expected_parameters = list(
            [[1, 2, 3, 78.67881492645, 128.2771922378],
             [1, 2, 7, 34.202963712735, 133.1339832262],
             [2, 1, 4, 50.00415699765, 110.81136453],
             [2, 1, 5, 50.00415699765, 110.81136453],
             [2, 1, 6, 50.00415699765, 110.81136453],
             [3, 2, 7, 34.202963712735, 133.1339832262],
             [4, 1, 5, 33.78875634641, 110.2468561538],
             [4, 1, 6, 33.78875634641, 110.2468561538],
             [5, 1, 6, 33.78875634641, 110.2468561538]])

        # Check resulting parameters
        for angle in topology.angles:
            w_angle = WritableAngle(angle)
            w_parameters = [attr[1] for attr in list(w_angle)]
            assert w_parameters in expected_parameters, \
                'Invalid writable angle parameters {}'.format(w_parameters)


class TestDihedrals(object):
    """
    It wraps all tests that involve dihedral parameters.
    """

    def test_OFF_parameters(self):
        """
        It checks the standard Open Force Field parameterization for
        dihedrals.
        """

        # Load molecule
        molecule = Molecule(smiles='C=CC(=O)O')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Check resulting parameters for proper torsions
        expected_propers = [(0, 1, 2, 3,
                             unit.Quantity(0.603518062312,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (0, 1, 2, 3,
                             unit.Quantity(0.5248455212365,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(0.0, unit.degree),
                             3, 1.0),
                            (0, 1, 2, 4,
                             unit.Quantity(0.9350453896311,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (1, 2, 4, 8,
                             unit.Quantity(2.529110648699,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (2, 1, 0, 5,
                             unit.Quantity(5.376019778605,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (2, 1, 0, 6,
                             unit.Quantity(5.376019778605,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (3, 2, 1, 7,
                             unit.Quantity(0.9350453896311,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (3, 2, 4, 8,
                             unit.Quantity(2.237928151469,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (3, 2, 4, 8,
                             unit.Quantity(1.23728649144,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(0.0, unit.degree),
                             1, 1.0),
                            (4, 2, 1, 7,
                             unit.Quantity(0.9350453896311,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (5, 0, 1, 7,
                             unit.Quantity(5.376019778605,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0),
                            (6, 0, 1, 7,
                             unit.Quantity(5.376019778605,
                                           unit.kilocalorie / unit.mole),
                             unit.Quantity(180.0, unit.degree),
                             2, 1.0)]

        assert len(expected_propers) == len(parameters['propers']), \
            'Unexpected number of proper torsions'

        for proper in parameters['propers']:
            proper_parameters = (proper['atom1_idx'], proper['atom2_idx'],
                                 proper['atom3_idx'], proper['atom4_idx'],
                                 proper['k'], proper['phase'],
                                 proper['periodicity'], proper['idivf'])

            assert proper_parameters in expected_propers, \
                'Unexpected proper torsion'

        # Check resulting parameters for improper torsions
        expected_impropers = [(0, 1, 2, 7,
                               unit.Quantity(1.1,
                                             unit.kilocalorie / unit.mole),
                               unit.Quantity(180.0, unit.degree),
                               2, 1),
                              (1, 0, 5, 6,
                               unit.Quantity(1.1,
                                             unit.kilocalorie / unit.mole),
                               unit.Quantity(180.0, unit.degree),
                               2, 1),
                              (1, 2, 3, 4,
                               unit.Quantity(10.5,
                                             unit.kilocalorie / unit.mole),
                               unit.Quantity(180.0, unit.degree),
                               2, 1)]

        assert len(expected_impropers) == \
            len(parameters['impropers']), \
            'Unexpected number of improper torsions'

        for improper in parameters['impropers']:
            improper_parameters = (improper['atom1_idx'],
                                   improper['atom2_idx'],
                                   improper['atom3_idx'],
                                   improper['atom4_idx'],
                                   improper['k'],
                                   improper['phase'],
                                   improper['periodicity'],
                                   improper['idivf'])

            assert improper_parameters in expected_impropers, \
                'Unexpected improper torsion'

    def test_Impact_writable_parameters(self):
        """
        It checks the Impact writable representation for dihedrals.
        """

        # Load benzene ring
        molecule = Molecule(smiles='CC=O')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Generate topology
        topology = Topology(molecule, parameters)

        expected_parameters = list(
            [[3, 2, 1, 4, 0.5107183999341, 1, 1, 0.0],
             [3, 2, 1, 5, 0.5107183999341, 1, 1, 0.0],
             [3, 2, 1, 6, 0.5107183999341, 1, 1, 0.0],
             [4, 1, 2, 7, 0.1493988458474, 1, 3, 0.0],
             [5, 1, 2, 7, 0.1493988458474, 1, 3, 0.0],
             [6, 1, 2, 7, 0.1493988458474, 1, 3, 0.0],
             [3, 2, 1, 4, -0.0275194074427, 1, 2, 0.0],
             [3, 2, 1, 5, -0.0275194074427, 1, 2, 0.0],
             [3, 2, 1, 6, -0.0275194074427, 1, 2, 0.0],
             [3, 2, 1, 4, -0.1057540923121, -1, 3, 0.0],
             [3, 2, 1, 5, -0.1057540923121, -1, 3, 0.0],
             [3, 2, 1, 6, -0.1057540923121, -1, 3, 0.0]])

        # Check resulting parameters
        for proper in topology.propers:
            w_proper = WritableProper(proper)
            w_parameters = [attr[1] for attr in list(w_proper)]
            assert w_parameters in expected_parameters, \
                'Invalid writable proper parameters {}'.format(w_parameters)

        expected_parameters = list([[1, 2, 3, 7, 1.1, -1, 2]])

        # Check resulting parameters
        for improper in topology.impropers:
            w_improper = WritableImproper(improper)
            w_parameters = [attr[1] for attr in list(w_improper)]
            assert w_parameters in expected_parameters, \
                'Invalid writable improper parameters {}'.format(w_parameters)

    def test_OFF_to_PELE_conversion(self):
        """
        It checks the difference between dihedral equations from PELE and
        Open Force Field. Their values should match throughout all the domain.
        """

        MAX_THRESHOLD = 1e-10

        for ligand_path in SET_OF_LIGAND_PATHS:
            ligand_path = get_data_file_path(ligand_path)

            # Load molecule
            molecule = Molecule(ligand_path)

            # Parameterize
            ff = OpenForceField(FORCEFIELD_NAME)
            parameters = ff.parameterize(molecule, charge_method='gasteiger')

            # Generate topology
            topology = Topology(molecule, parameters)

            x = unit.Quantity(np.arange(0, np.pi, 0.1), unit=unit.radians)

            for PELE_proper, OFF_proper in zip(topology.propers,
                                               topology._OFF_propers):
                PELE_y = apply_PELE_dihedral_equation(PELE_proper, x)
                OFF_y = apply_OFF_dihedral_equation(OFF_proper, x)

                y_diff = PELE_y - OFF_y

                assert np.linalg.norm(y_diff) < MAX_THRESHOLD

    def test_excluded_dihedrals_handler(self):
        """
        It checks the exclusion of dihedrals when generating 1-4 lists
        for PELE. This process is done internally in PELE but we
        need to mark previously those dihedrals to exclude from 1-4
        lists.
        """

        # Load molecule
        molecule = Molecule()

        # Parameterize
        parameters = BaseParameterWrapper()

        # Generate topology
        topology = Topology(molecule, parameters)

        # Add several dummy propers
        topology.add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=3,
                   periodicity=1, prefactor=1, index=0,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        topology.add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=3,
                   periodicity=2, prefactor=1, index=1,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        topology.add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=5, atom4_idx=3,
                   periodicity=2, prefactor=1, index=2,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        topology.add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=4,
                   periodicity=1, prefactor=1, index=3,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        topology.add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=5,
                   periodicity=1, prefactor=1, index=4,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        topology.add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=6,
                   periodicity=1, prefactor=1, index=5,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))

        # Add one bond
        topology.add_bond(
            Bond(atom1_idx=0, atom2_idx=4,
                 spring_constant=unit.Quantity(1.0, unit.kilocalorie
                                               / (unit.angstrom ** 2
                                                  * unit.mole)),
                 eq_dist=unit.Quantity(1.0, unit.angstrom)))

        # Add one angle
        topology.add_angle(
            Angle(atom1_idx=0, atom2_idx=1, atom3_idx=5,
                  spring_constant=unit.Quantity(1.0, unit.kilocalorie
                                                / (unit.radian ** 2
                                                   * unit.mole)),
                  eq_angle=unit.Quantity(1.0, unit.degrees)))

        topology._handle_excluded_propers()

        # Assertions
        for proper in topology.propers:
            if proper.index == 0:
                assert proper.exclude is False, \
                    'Proper 0 must be INCLUDED in 1-4 list'
            elif proper.index == 1:
                assert proper.exclude is False, \
                    'Proper 1 must be INCLUDED in 1-4 list'
            elif proper.index == 2:
                assert proper.exclude is True, \
                    'Proper 2 must be EXCLUDED in 1-4 list'
            elif proper.index == 3:
                assert proper.exclude is True, \
                    'Proper 3 must be EXCLUDED in 1-4 list'
            elif proper.index == 4:
                assert proper.exclude is True, \
                    'Proper 4 must be EXCLUDED in 1-4 list'
            elif proper.index == 5:
                assert proper.exclude is False, \
                    'Proper 5 must be INCLUDED in 1-4 list'

    def test_nonstandard_dihedrals_writable_parameters(self):
        """
        It checks the writable representation of non standard dihedrals.
        """
        # Load molecule
        molecule = Molecule(smiles='c1c(c(n(n1)S(=O)(=O)C))O')

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Generate topology
        topology = Topology(molecule, parameters)

        expected_parameters = [[1, 2, -3, 4, 5.376019778605, -1, 2, 0.0],
                               [1, 2, 3, 12, 5.376019778605, -1, 2, 0.0],
                               [1, 2, 10, 16, 0.8722932201352, -1, 2, 0.0],
                               [1, 5, -4, 3, -0.1847566338874, 1, 3, 0.0],
                               [1, 5, 4, 6, -0.1847566338874, 1, 3, 0.0],
                               [2, 1, -5, 4, 6.867249875207, -1, 2, 0.0],
                               [2, 3, -4, 5, 1.131350927723, -1, 2, 0.0],
                               [2, 3, 4, 6, 1.131350927723, -1, 2, 0.0],
                               [3, 2, -1, 5, 0.9350453896311, -1, 2, 0.0],
                               [3, 2, 1, 11, 0.9350453896311, -1, 2, 0.0],
                               [3, 2, 10, 16, 0.8722932201352, -1, 2, 0.0],
                               [3, 4, 6, 7, 0.5146229154488, 1, 1, 0.0],
                               [3, 4, 6, 8, 0.5146229154488, 1, 1, 0.0],
                               [3, 4, 6, 9, -0.5050335923881, 1, 3, 90.0],
                               [4, 3, 2, 10, 5.376019778605, -1, 2, 0.0],
                               [4, 5, 1, 11, 6.867249875207, -1, 2, 0.0],
                               [4, 6, 9, 13, 0.1055132590818, 1, 3, 0.0],
                               [4, 6, 9, 14, 0.1055132590818, 1, 3, 0.0],
                               [4, 6, 9, 15, 0.1055132590818, 1, 3, 0.0],
                               [5, 1, 2, 10, 0.9350453896311, -1, 2, 0.0],
                               [5, 4, 3, 12, 1.131350927723, -1, 2, 0.0],
                               [5, 4, 6, 7, 0.082636589537, 1, 1, 0.0],
                               [5, 4, 6, 8, 0.082636589537, 1, 1, 0.0],
                               [5, 4, 6, 9, 1.914473924133, 1, 1, 0.0],
                               [6, 4, 3, 12, 1.131350927723, -1, 2, 0.0],
                               [7, 6, 9, 13, 0.1055132590818, 1, 3, 0.0],
                               [7, 6, 9, 14, 0.1055132590818, 1, 3, 0.0],
                               [7, 6, 9, 15, 0.1055132590818, 1, 3, 0.0],
                               [8, 6, 9, 13, 0.1055132590818, 1, 3, 0.0],
                               [8, 6, 9, 14, 0.1055132590818, 1, 3, 0.0],
                               [8, 6, 9, 15, 0.1055132590818, 1, 3, 0.0],
                               [10, 2, 1, 11, 0.9350453896311, -1, 2, 0.0],
                               [10, 2, 3, 12, 5.376019778605, -1, 2, 0.0],
                               [1, 5, -4, 3, 1.471695489046, -1, 2, 0.0],
                               [1, 5, 4, 6, 1.471695489046, -1, 2, 0.0],
                               [3, 4, 6, 9, 0.3649469393198, 1, 2, 0.0]]

        # Check resulting parameters
        for proper in topology.propers:
            w_proper = WritableProper(proper)
            w_parameters = [attr[1] for attr in list(w_proper)]
            assert w_parameters in expected_parameters, \
                'Invalid writable proper parameters {}'.format(w_parameters)

        expected_parameters = [[1, 2, 3, 10, 1.1, -1, 2],
                               [2, 1, 5, 11, 1.1, -1, 2],
                               [2, 3, 4, 12, 1.1, -1, 2],
                               [3, 4, 5, 6, 10.5, -1, 2]]

        # Check resulting parameters
        for improper in topology.impropers:
            w_improper = WritableImproper(improper)
            w_parameters = [attr[1] for attr in list(w_improper)]
            assert w_parameters in expected_parameters, \
                'Invalid writable improper parameters {}'.format(w_parameters)

    def test_OPLS_dummy_propers(self):
        """
        It checks that the number of propers that are obtained from OPLS
        is correct, ensuring that dummy propers (propers with a null
        force constant whose definitions are only used in the 1-4 lists
        of PELE) are present in the parameterized proper list.
        """
        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper

        # Load molecule
        molecule = Molecule(
            get_data_file_path('ligands/octafluorocyclobutane.pdb'))

        # Set force field and obtain parameters
        ffld_file = get_data_file_path('tests/CO1_ffld_output.txt')
        with open(ffld_file) as f:
            ffld_output = f.read()

        parameters = OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                               ffld_output)

        # Generate topology
        topology = Topology(molecule, parameters)

        # Validate number of propers
        assert len(topology.propers) == 100, \
            'Unexpected number of proper torsions'


class TestCharges(object):
    """
    It wraps all tests that involve charge parameters.
    """
    LIGAND_PATH = 'ligands/oleic_acid.pdb'

    def test_am1bcc_method(self):
        """It tests the am1bcc method"""
        ligand_path = get_data_file_path(self.LIGAND_PATH)

        # Load molecule
        molecule = Molecule(ligand_path)

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='am1bcc')

        # Check charges
        check_CHO_charges(parameters)

    def test_gasteiger_method(self):
        """It tests the gasteiger method"""
        ligand_path = get_data_file_path(self.LIGAND_PATH)

        # Load molecule
        molecule = Molecule(ligand_path)

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='gasteiger')

        # Check charges
        check_CHO_charges(parameters)

    def test_OPLS_method(self):
        """It tests the OPLS method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)

        # Load molecule
        molecule = Molecule(ligand_path)

        expected_charges = [-0.22, 0.7, -0.12, -0.8, -0.8, -0.12, -0.12,
                            -0.12, -0.12, -0.12, -0.115, -0.115, -0.12,
                            -0.12, -0.12, -0.12, -0.12, -0.12, -0.12,
                            -0.18, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
                            0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
                            0.06, 0.115, 0.115, 0.06, 0.06, 0.06, 0.06,
                            0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
                            0.06, 0.06, 0.06, 0.06, 0.06, 0.06]

        # Workaround to avoid the use of the Schrodinger Toolkit
        from .utils import parameterize_opls2005

        # Load OPLS2005 force field and locate ffld_file
        oplsff = OPLS2005ForceField()
        ffld_file = get_data_file_path('tests/OLC_ffld_output.txt')

        # Parameterize
        parameters = parameterize_opls2005(oplsff, molecule, ffld_file)

        partial_charges = parameters['charges']

        assert len(partial_charges) == len(expected_charges), \
            'Size of Molecule\'s partial charges is expected to match ' \
            + 'with size of reference charges list'

        for partial_charges, expected_charge in zip(partial_charges,
                                                    expected_charges):
            assert partial_charges == unit.Quantity(expected_charge,
                                                    unit.elementary_charge), \
                'Unexpected partial charge {}'.format(partial_charges)

    def test_dummy_method(self):
        """It tests the dummy charge calculator."""

        ligand_path = get_data_file_path(self.LIGAND_PATH)

        # Load molecule
        molecule = Molecule(ligand_path)

        expected_charges = [0.0, ] * 53

        # Parameterize
        ff = OpenForceField(FORCEFIELD_NAME)
        parameters = ff.parameterize(molecule, charge_method='dummy')

        partial_charges = parameters['charges']

        assert len(partial_charges) == len(expected_charges), \
            'Size of Molecule\'s partial charges is expected to match ' \
            + 'with size of reference charges list'

        for partial_charges, expected_charge in zip(partial_charges,
                                                    expected_charges):
            assert partial_charges == unit.Quantity(expected_charge,
                                                    unit.elementary_charge), \
                'Unexpected partial charge {}'.format(partial_charges)


class TestSolventParameters(object):
    """
    It holds tests to validate the assignment of solvent parameters.
    """

    def test_add_SGBNP_solvent_parameters(self):
        """
        It tests the function that adds the SGBNP solvent parameters to
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

        OPLS2005ParameterWrapper._add_SGBNP_solvent_parameters(params1)
        OPLS2005ParameterWrapper._add_SGBNP_solvent_parameters(params2)
        OPLS2005ParameterWrapper._add_SGBNP_solvent_parameters(params3)

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

    def test_GBSA_params_by_type(self):
        """
        It tests the params_by_type dictionary in the assignment of the
        GBSA parameters.
        """

        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper

        # 1st test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' C1 ']
        OPLS_params['atom_types'] = [' CT']
        degree_by_name = dict(((' C1 ', 4), ))
        parent_by_name = dict(((' C1 ', None), ))
        element_by_name = dict(((' C1 ', 'C'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.72, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.9], \
            'Unexpected GBSA scale'

        # 2nd test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' C1 ']
        OPLS_params['atom_types'] = [' CW']
        degree_by_name = dict(((' C1 ', 3), ))
        parent_by_name = dict(((' C1 ', None), ))
        element_by_name = dict(((' C1 ', 'C'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.72, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.875], \
            'Unexpected GBSA scale'

        # 3rd test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' O1 ']
        OPLS_params['atom_types'] = [' O ']
        degree_by_name = dict(((' O1 ', 1), ))
        parent_by_name = dict(((' O1 ', None), ))
        element_by_name = dict(((' O1 ', 'O'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.48], \
            'Unexpected GBSA scale'

        # 4th test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' H10']
        OPLS_params['atom_types'] = [' H2']
        degree_by_name = dict(((' H10', 1), ))
        parent_by_name = dict(((' H10', 'C'), ))
        element_by_name = dict(((' H10', 'H'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.25], \
            'Unexpected GBSA scale'

    def test_GBSA_params_by_element(self):
        """
        It tests the scale_by_element and radius_by_element dictionaries
        in the assignment of the GBSA parameters.
        """

        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper

        # 1st test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' C1 ']
        OPLS_params['atom_types'] = [' C?!']
        degree_by_name = dict(((' C1 ', 3), ))
        parent_by_name = dict(((' C1 ', None), ))
        element_by_name = dict(((' C1 ', 'C'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.72, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.875], \
            'Unexpected GBSA scale'

        # 2nd test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' O1 ']
        OPLS_params['atom_types'] = ['O?!']
        degree_by_name = dict(((' O1 ', 2), ))
        parent_by_name = dict(((' O1 ', None), ))
        element_by_name = dict(((' O1 ', 'O'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.535], \
            'Unexpected GBSA scale'

        # 3rd test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' O1 ']
        OPLS_params['atom_types'] = ['O?!']
        degree_by_name = dict(((' O1 ', 1), ))
        parent_by_name = dict(((' O1 ', None), ))
        element_by_name = dict(((' O1 ', 'O'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.48], \
            'Unexpected GBSA scale'

        # 4th test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' H10']
        OPLS_params['atom_types'] = ['H?!']
        degree_by_name = dict(((' H10', 1), ))
        parent_by_name = dict(((' H10', 'C'), ))
        element_by_name = dict(((' H10', 'H'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.25], \
            'Unexpected GBSA scale'

        # 5th test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' H10']
        OPLS_params['atom_types'] = ['H?!']
        degree_by_name = dict(((' H10', 1), ))
        parent_by_name = dict(((' H10', 'N'), ))
        element_by_name = dict(((' H10', 'H'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.15], \
            'Unexpected GBSA scale'

        # 6th test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' H10']
        OPLS_params['atom_types'] = ['H?!']
        degree_by_name = dict(((' H10', 1), ))
        parent_by_name = dict(((' H10', 'O'), ))
        element_by_name = dict(((' H10', 'H'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.85, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.05], \
            'Unexpected GBSA scale'

        # 7th test
        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' C1 ']
        OPLS_params['atom_types'] = [' C?!']
        degree_by_name = dict(((' C1 ', 2), ))
        parent_by_name = dict(((' C1 ', None), ))
        element_by_name = dict(((' C1 ', 'C'), ))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(OPLS_params,
                                                              degree_by_name,
                                                              parent_by_name,
                                                              element_by_name)

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.72, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [1.825], \
            'Unexpected GBSA scale'

    def test_GBSA_default_params(self):
        """
        It tests the default parameters of the GBSA implementation.
        """
        import io
        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper
        from peleffy.utils import Logger

        OPLS_params = OPLS2005ParameterWrapper()

        # Create mock molecule containing just one customized atom
        OPLS_params['atom_names'] = [' C1 ']
        OPLS_params['atom_types'] = [' C?!']
        degree_by_name = dict(((' C1 ', 2), ))
        parent_by_name = dict(((' C1 ', None), ))
        element_by_name = dict(((' C1 ', '?'), ))

        import logging

        # Force a hard reset of logging library and the logger it manages
        from importlib import reload
        logging.shutdown()
        reload(logging)

        # Initiate logger
        log = Logger()

        # Try the default level (INFO)
        # Catch logger messages to string buffer
        with io.StringIO() as buf:
            # Add custom handler to logger
            log_handler = logging.StreamHandler(buf)
            log._logger.handlers = list()
            log._logger.addHandler(log_handler)

            OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(
                OPLS_params, degree_by_name, parent_by_name, element_by_name)

            # Get string from buffer
            output = buf.getvalue()

            assert output == 'Warning: OBC parameters for ' \
                + 'C1 C?! NOT found in the template database. ' \
                + 'Using default parameters\n'

        assert OPLS_params['GBSA_radii'] == \
            [unit.Quantity(0.80, unit.angstrom)], \
            'Unexpected GBSA radii'
        assert OPLS_params['GBSA_scales'] == [2.0], \
            'Unexpected GBSA scale'

    def test_add_GBSA_solvent_parameters(self):
        """
        It tests the function that adds the GBSA solvent parameters to
        the OPLSParameters collection.
        """

        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper
        from peleffy.topology import Molecule
        from peleffy.utils import get_data_file_path

        def generate_and_compare_GBSA_solvent_parameters(pdb_file,
                                                         ffld_file,
                                                         reference_txt):
            """
            Given a ligand, it tests that the HTC radii and the scale factor
            for heteroatoms computed with the function
            _add_GBSA_solvent_parameters correspond to the ones
            obtained with the original solventOBCParamsGenerator.py script.

            Parameters
            ----------
            pdb_file : str
                The path to the PDB of the ligand to test
            ffld_file : str
                The path to the ffld_server's output file
            reference_txt : str
                The path to reference TXT file obtained
                from solventOBCParamsGenerator.py
            """

            # Load the molecule
            molecule = Molecule(pdb_file)

            # Set force field and obtain parameters
            with open(ffld_file) as f:
                ffld_output = f.read()

            # Initializate wrapper for OPLS2005 parameters
            wrapper_opls = OPLS2005ParameterWrapper()
            parameters = wrapper_opls.from_ffld_output(molecule, ffld_output)

            # Get the radi and scale parameters for each atom name
            atom_names = [param.replace('_', '')
                          for param in parameters['atom_names']]
            d_radii = dict(zip(atom_names, parameters['GBSA_radii']))
            d_scale = dict(zip(atom_names, parameters['GBSA_scales']))

            # Load the reference file obtained from solventOBCParamsGenerator.py
            data = open(reference_txt, 'r').readlines()

            # Check that the parameters correspond
            for line in data:
                params = line.split()
                assert float(params[3]) == d_scale.get(params[1]), \
                    'Unexpected GBSA Overlap factor'
                assert unit.Quantity(float(params[4]), unit.angstrom) == \
                    d_radii.get(params[1]), \
                    'Unexpected GBSA radius'

        # For malonate
        ffld_file = get_data_file_path('tests/MAL_ffld_output.txt')
        pdb_file = get_data_file_path('ligands/malonate.pdb')
        reference_txt = get_data_file_path('tests/malz_OBCParams.txt')
        generate_and_compare_GBSA_solvent_parameters(pdb_file,
                                                     ffld_file,
                                                     reference_txt)

        # For methane
        ffld_file = get_data_file_path('tests/MET_ffld_output.txt')
        pdb_file = get_data_file_path('ligands/methane.pdb')
        reference_txt = get_data_file_path('tests/metz_OBCParams.txt')
        generate_and_compare_GBSA_solvent_parameters(pdb_file,
                                                     ffld_file,
                                                     reference_txt)

        # For ethylene
        ffld_file = get_data_file_path('tests/ETL_ffld_output.txt')
        pdb_file = get_data_file_path('ligands/ethylene.pdb')
        reference_txt = get_data_file_path('tests/etlz_OBCParams.txt')
        generate_and_compare_GBSA_solvent_parameters(pdb_file,
                                                     ffld_file,
                                                     reference_txt)
