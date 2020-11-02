"""
This module contains the tests to check peleffy's parameters.
"""

import pytest

from simtk import unit
import numpy as np

from peleffy.utils import get_data_file_path
from .utils import (SET_OF_LIGAND_PATHS, apply_PELE_dihedral_equation,
                    apply_OFF_dihedral_equation, check_CHO_charges_in_molecule)
from peleffy.topology import Molecule, Bond, Angle, Proper
from peleffy.template.impact import (WritableBond, WritableAngle,
                                     WritableProper, WritableImproper)


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'


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
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

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

        for bond in molecule.parameters['bonds']:
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
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

        expected_parameters = list([[1, 2, 332.5750972667, 1.523640340452],
                                    [1, 4, 376.8940758588, 1.094223427522],
                                    [1, 5, 376.8940758588, 1.094223427522],
                                    [1, 6, 376.8940758588, 1.094223427522],
                                    [2, 3, 608.3286693405, 1.225108345696],
                                    [2, 7, 404.20804685, 1.085503378387]])

        # Check resulting parameters
        for bond in molecule.bonds:
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
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

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

        for angle in molecule.parameters['angles']:
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
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

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
        for angle in molecule.angles:
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

        MAX_THRESHOLD = 1e-3

        # Load molecule
        molecule = Molecule(smiles='C=CC(=O)O')

        # Parameterize
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

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

        assert len(expected_propers) == len(molecule.parameters['propers']), \
            'Unexpected number of proper torsions'

        for proper in molecule.parameters['propers']:
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
            len(molecule.parameters['impropers']), \
            'Unexpected number of improper torsions'

        for improper in molecule.parameters['impropers']:
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
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

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
        for proper in molecule.propers:
            w_proper = WritableProper(proper)
            w_parameters = [attr[1] for attr in list(w_proper)]
            assert w_parameters in expected_parameters, \
                'Invalid writable proper parameters {}'.format(w_parameters)

        expected_parameters = list([[1, 2, 3, 7, 1.1, -1, 2]])

        # Check resulting parameters
        for improper in molecule.impropers:
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
            molecule = Molecule(ligand_path)
            molecule.parameterize(FORCEFIELD_NAME, charge_method='gasteiger')

            x = unit.Quantity(np.arange(0, np.pi, 0.1), unit=unit.radians)

            for PELE_proper, OFF_proper in zip(molecule.propers,
                                               molecule._OFF_propers):
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
        molecule = Molecule()

        # Add several dummy propers
        molecule._add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=3,
                   periodicity=1, prefactor=1, index=0,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        molecule._add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=3,
                   periodicity=2, prefactor=1, index=1,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        molecule._add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=5, atom4_idx=3,
                   periodicity=2, prefactor=1, index=2,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        molecule._add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=4,
                   periodicity=1, prefactor=1, index=3,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        molecule._add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=5,
                   periodicity=1, prefactor=1, index=4,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))
        molecule._add_proper(
            Proper(atom1_idx=0, atom2_idx=1, atom3_idx=2, atom4_idx=6,
                   periodicity=1, prefactor=1, index=5,
                   constant=unit.Quantity(1.0, unit.kilocalorie / unit.mole)))

        # Add one bond
        molecule._add_bond(
            Bond(atom1_idx=0, atom2_idx=4,
                 spring_constant=unit.Quantity(1.0, unit.kilocalorie
                                               / (unit.angstrom ** 2
                                                  * unit.mole)),
                 eq_dist=unit.Quantity(1.0, unit.angstrom)))

        # Add one angle
        molecule._add_angle(
            Angle(atom1_idx=0, atom2_idx=1, atom3_idx=5,
                  spring_constant=unit.Quantity(1.0, unit.kilocalorie
                                                / (unit.radian ** 2
                                                   * unit.mole)),
                  eq_angle=unit.Quantity(1.0, unit.degrees)))

        molecule._handle_excluded_propers()

        # Assertions
        for proper in molecule.propers:
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
        molecule = Molecule(smiles='c1c(c(n(n1)S(=O)(=O)C))O')

        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charge_method='gasteiger')

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
        for proper in molecule.propers:
            w_proper = WritableProper(proper)
            w_parameters = [attr[1] for attr in list(w_proper)]
            assert w_parameters in expected_parameters, \
                'Invalid writable proper parameters {}'.format(w_parameters)

        expected_parameters = [[1, 2, 3, 10, 1.1, -1, 2],
                               [2, 1, 5, 11, 1.1, -1, 2],
                               [2, 3, 4, 12, 1.1, -1, 2],
                               [3, 4, 5, 6, 10.5, -1, 2]]

        # Check resulting parameters
        for improper in molecule.impropers:
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
        from peleffy.forcefield import OPLS2005ForceField
        from peleffy.forcefield import OPLS2005ParameterWrapper

        # Load molecule
        molecule = Molecule(get_data_file_path('ligands/CO1.pdb'))
        oplsff = OPLS2005ForceField('OPLS2005')

        # Set force field and obtain parameters
        ffld_file = get_data_file_path('tests/CO1_ffld_output.txt')
        with open(ffld_file) as f:
            ffld_output = f.read()

        parameters = OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                               ffld_output)
        oplsff._parameters = parameters

        molecule.set_forcefield(oplsff)
        molecule.parameterize(charge_method='gasteiger')

        # Validate number of propers
        assert len(molecule.propers) == 100, \
            'Unexpected number of proper torsions'


class TestCharges(object):
    """
    It wraps all tests that involve charge parameters.
    """
    LIGAND_PATH = 'ligands/OLC.pdb'

    def test_am1bcc_method(self):
        """It tests the am1bcc method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)
        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME, charge_method='am1bcc')
        check_CHO_charges_in_molecule(molecule)

    def test_gasteiger_method(self):
        """It tests the gasteiger method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)
        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME, charge_method='gasteiger')
        check_CHO_charges_in_molecule(molecule)

    def test_OPLS_method(self):
        """It tests the OPLS method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)
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
        from peleffy.forcefield import OPLS2005ParameterWrapper

        ffld_file = get_data_file_path('tests/OLC_ffld_output.txt')
        with open(ffld_file) as f:
            ffld_output = f.read()
        molecule._parameters = \
            OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                      ffld_output)

        # Run charge calculator
        from peleffy.charge import OPLSChargeCalculator
        charge_calculator = OPLSChargeCalculator(molecule)
        partial_charges = charge_calculator.get_partial_charges()

        assert len(partial_charges) == len(expected_charges), \
            'Size of Molecule\'s partial charges is expected to match ' \
            + 'with size of reference charges list'

        for partial_charges, expected_charge in zip(partial_charges,
                                                    expected_charges):
            assert partial_charges == unit.Quantity(expected_charge,
                                                    unit.elementary_charge), \
                'Unexpected partial charge {}'.format(partial_charges)
