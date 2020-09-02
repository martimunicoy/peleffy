"""
This module contains the tests to check offpele's parameters.
"""

import pytest

from simtk import unit
import numpy as np

from offpele.utils import get_data_file_path
from offpele.utils.toolkits import SchrodingerToolkitWrapper
from .utils import (SET_OF_LIGAND_PATHS, apply_PELE_dihedral_equation,
                    apply_OFF_dihedral_equation, check_CHO_charges_in_molecule)
from offpele.topology import Molecule, Bond, Angle, Proper
from offpele.template.impact import (WritableBond, WritableAngle,
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
                              charges_method='gasteiger')

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

        for indexes, properties in dict(molecule.parameters['Bonds']).items():
            expected_length = unit.Quantity(expected_lengths[indexes],
                                            unit.angstrom)
            expected_k = unit.Quantity(expected_ks[indexes],
                                       unit.kilocalorie
                                       / (unit.angstrom ** 2 * unit.mole))

            assert properties.length - expected_length \
                < unit.Quantity(MAX_THRESHOLD, unit.angstrom), \
                'Invalid length for bond {} {}'.format(indexes, properties)
            assert properties.k - expected_k \
                < unit.Quantity(MAX_THRESHOLD, unit.kilocalorie
                                / (unit.angstrom ** 2 * unit.mole)), \
                'Invalid k for bond {} {}'.format(indexes, properties)

    def test_Impact_writable_parameters(self):
        """
        It checks the Impact writable representation for bonds.
        """

        # Load benzene ring
        molecule = Molecule(smiles='CC=O')

        # Parameterize
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charges_method='gasteiger')

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
                              charges_method='gasteiger')

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

        for indexes, properties in dict(molecule.parameters['Angles']).items():
            expected_angle = unit.Quantity(expected_angles[indexes],
                                           unit.degree)
            expected_k = unit.Quantity(expected_ks[indexes],
                                       unit.kilocalorie
                                       / (unit.radian ** 2 * unit.mole))

            assert properties.angle - expected_angle \
                < unit.Quantity(MAX_THRESHOLD, unit.degree), \
                'Invalid length for angle {} {}'.format(indexes, properties)
            assert properties.k - expected_k \
                < unit.Quantity(MAX_THRESHOLD, unit.kilocalorie
                                / (unit.radian ** 2 * unit.mole)), \
                'Invalid k for angle {} {}'.format(indexes, properties)

    def test_Impact_writable_parameters(self):
        """
        It checks the Impact writable representation for angles.
        """

        # Load benzene ring
        molecule = Molecule(smiles='CC=O')

        # Parameterize
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charges_method='gasteiger')

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
                              charges_method='gasteiger')

        # Check resulting parameters for proper torsions
        expected_ks = {(0, 1, 2, 3): [0.603518062312, 0.5248455212365],
                       (0, 1, 2, 4): [0.9350453896311],
                       (1, 2, 4, 8): [2.529110648699],
                       (2, 1, 0, 5): [5.376019778605],
                       (2, 1, 0, 6): [5.376019778605],
                       (3, 2, 1, 7): [0.9350453896311],
                       (3, 2, 4, 8): [2.237928151469, 1.23728649144],
                       (4, 2, 1, 7): [0.9350453896311],
                       (5, 0, 1, 7): [5.376019778605],
                       (6, 0, 1, 7): [5.376019778605]}

        expected_phases = {(0, 1, 2, 3): [180.0, 0.0],
                           (0, 1, 2, 4): [180.0],
                           (1, 2, 4, 8): [180.0],
                           (2, 1, 0, 5): [180.0],
                           (2, 1, 0, 6): [180.0],
                           (3, 2, 1, 7): [180.0],
                           (3, 2, 4, 8): [180.0, 0.0],
                           (4, 2, 1, 7): [180.0],
                           (5, 0, 1, 7): [180.0],
                           (6, 0, 1, 7): [180.0]}

        expected_periodicities = {(0, 1, 2, 3): [2, 3],
                                  (0, 1, 2, 4): [2],
                                  (1, 2, 4, 8): [2],
                                  (2, 1, 0, 5): [2],
                                  (2, 1, 0, 6): [2],
                                  (3, 2, 1, 7): [2],
                                  (3, 2, 4, 8): [2, 1],
                                  (4, 2, 1, 7): [2],
                                  (5, 0, 1, 7): [2],
                                  (6, 0, 1, 7): [2]}

        expected_idivfs = {(0, 1, 2, 3): [1.0, 1.0],
                           (0, 1, 2, 4): [1.0],
                           (1, 2, 4, 8): [1.0],
                           (2, 1, 0, 5): [1.0],
                           (2, 1, 0, 6): [1.0],
                           (3, 2, 1, 7): [1.0],
                           (3, 2, 4, 8): [1.0, 1.0],
                           (4, 2, 1, 7): [1.0],
                           (5, 0, 1, 7): [1.0],
                           (6, 0, 1, 7): [1.0]}

        for indexes, properties in dict(
                molecule.parameters['ProperTorsions']).items():
            for i, (k, phase, periodicity, idivf) in enumerate(
                    zip(properties.k, properties.phase,
                        properties.periodicity, properties.idivf)):
                expected_k = unit.Quantity(expected_ks[indexes][i],
                                           unit.kilocalorie / unit.mole)
                expected_phase = unit.Quantity(expected_phases[indexes][i],
                                               unit.degree)
                expected_periodicity = expected_periodicities[indexes][i]
                expected_idivf = expected_idivfs[indexes][i]

                assert k - expected_k < \
                    unit.Quantity(MAX_THRESHOLD,
                                  unit.kilocalorie / unit.mole), \
                    'Invalid k for proper torsion ' \
                    + '{} {}'.format(indexes, properties)

                assert phase - expected_phase < \
                    unit.Quantity(MAX_THRESHOLD, unit.degree), \
                    'Invalid phase for proper torsion ' \
                    + '{} {}'.format(indexes, properties)

                assert periodicity - expected_periodicity < MAX_THRESHOLD, \
                    'Invalid periodicity for proper torsion ' \
                    + '{} {}'.format(indexes, properties)

                assert idivf - expected_idivf < MAX_THRESHOLD, \
                    'Invalid idivf for proper torsion ' \
                    + '{} {}'.format(indexes, properties)

        # Check resulting parameters for improper torsions
        expected_ks = {(0, 1, 2, 7): [1.1],
                       (1, 0, 5, 6): [1.1],
                       (1, 2, 3, 4): [10.5]}

        expected_phases = {(0, 1, 2, 7): [180.0],
                           (1, 0, 5, 6): [180.0],
                           (1, 2, 3, 4): [180.0]}

        expected_periodicities = {(0, 1, 2, 7): [2],
                                  (1, 0, 5, 6): [2],
                                  (1, 2, 3, 4): [2]}

        for indexes, properties in dict(
                molecule.parameters['ImproperTorsions']).items():
            for i, (k, phase, periodicity) in enumerate(
                    zip(properties.k, properties.phase,
                        properties.periodicity)):
                expected_k = unit.Quantity(expected_ks[indexes][i],
                                           unit.kilocalorie / unit.mole)
                expected_phase = unit.Quantity(expected_phases[indexes][i],
                                               unit.degree)
                expected_periodicity = expected_periodicities[indexes][i]

                assert k - expected_k < \
                    unit.Quantity(MAX_THRESHOLD,
                                  unit.kilocalorie / unit.mole), \
                    'Invalid k for improper torsion ' \
                    + '{} {}'.format(indexes, properties)

                assert phase - expected_phase < \
                    unit.Quantity(MAX_THRESHOLD, unit.degree), \
                    'Invalid phase for improper torsion ' \
                    + '{} {}'.format(indexes, properties)

                assert periodicity - expected_periodicity < MAX_THRESHOLD, \
                    'Invalid periodicity for improper torsion ' \
                    + '{} {}'.format(indexes, properties)

            assert properties.idivf is None, \
                'Invalid idivf for improper torsion ' \
                + '{} {}'.format(indexes, properties)

    def test_Impact_writable_parameters(self):
        """
        It checks the Impact writable representation for dihedrals.
        """

        # Load benzene ring
        molecule = Molecule(smiles='CC=O')

        # Parameterize
        molecule.parameterize('openff_unconstrained-1.2.0.offxml',
                              charges_method='gasteiger')

        expected_parameters = list(
            [[3, 2, 1, 4, 0.5107183999341, 1, 1],
             [3, 2, 1, 5, 0.5107183999341, 1, 1],
             [3, 2, 1, 6, 0.5107183999341, 1, 1],
             [4, 1, 2, 7, 0.1493988458474, 1, 3],
             [5, 1, 2, 7, 0.1493988458474, 1, 3],
             [6, 1, 2, 7, 0.1493988458474, 1, 3],
             [3, 2, 1, 4, -0.0275194074427, 1, 2],
             [3, 2, 1, 5, -0.0275194074427, 1, 2],
             [3, 2, 1, 6, -0.0275194074427, 1, 2],
             [3, 2, 1, 4, -0.1057540923121, -1, 3],
             [3, 2, 1, 5, -0.1057540923121, -1, 3],
             [3, 2, 1, 6, -0.1057540923121, -1, 3]])

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
            molecule.parameterize(FORCEFIELD_NAME, charges_method='gasteiger')

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


class TestCharges(object):
    """
    It wraps all tests that involve charge parameters.
    """
    LIGAND_PATH = 'ligands/OLC.pdb'

    def test_am1bcc_method(self):
        """It tests the am1bcc method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)
        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME, charges_method='am1bcc')
        check_CHO_charges_in_molecule(molecule)

    def test_gasteiger_method(self):
        """It tests the gasteiger method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)
        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME, charges_method='gasteiger')
        check_CHO_charges_in_molecule(molecule)

    def test_OPLS_method(self):
        """It tests the OPLS method"""

        ligand_path = get_data_file_path(self.LIGAND_PATH)
        molecule = Molecule(ligand_path)

        # To avoid the use of Schrodinger Toolkit
        charges = [-0.22, 0.7, -0.12, -0.8, -0.8, -0.12, -0.12, -0.12,
                   -0.12, -0.12, -0.115, -0.115, -0.12, -0.12, -0.12,
                   -0.12, -0.12, -0.12, -0.12, -0.18, 0.06, 0.06, 0.06,
                   0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
                   0.06, 0.06, 0.115, 0.115, 0.06, 0.06, 0.06, 0.06, 0.06,
                   0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06,
                   0.06, 0.06, 0.06]

        molecule._OPLS_parameters = SchrodingerToolkitWrapper.OPLSParameters(
            {'charges': [unit.Quantity(charge, unit.elementary_charge)
                         for charge in charges]})

        molecule.parameterize(FORCEFIELD_NAME, charges_method='OPLS')

        assert len(molecule.off_molecule.partial_charges) == len(charges), \
            'Size of Molecule\'s partial charges is expected to match ' \
            + 'with size of reference charges list'

        for charge, expected_charge in zip(
                molecule.off_molecule.partial_charges, charges):
            assert charge == unit.Quantity(expected_charge,
                                           unit.elementary_charge), \
                'Unexpected charge {}'.format(charge)
