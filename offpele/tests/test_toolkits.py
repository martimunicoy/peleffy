"""
This module contains the tests to check the offpele's toolkits.
"""

import pytest

import os

from offpele.topology import Molecule
from offpele.utils.toolkits import (SchrodingerToolkitWrapper,
                                    ToolkitUnavailableException)
from offpele.template.impact import (WritableAtom, WritableBond,
                                     WritableAngle)

from simtk import unit


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'
METHANE_OPLS_PARAMETERS = SchrodingerToolkitWrapper.OPLSParameters({
    'atom_names': [' C1 ', ' H2 ', ' H3 ', ' H4 ', ' H5 '],
    'atom_types': ['CT', 'HC', 'HC', 'HC', 'HC'],
    'charges': [unit.Quantity(-0.24, unit.elementary_charge),
                unit.Quantity(0.06, unit.elementary_charge),
                unit.Quantity(0.06, unit.elementary_charge),
                unit.Quantity(0.06, unit.elementary_charge),
                unit.Quantity(0.06, unit.elementary_charge)],
    'sigmas': [unit.Quantity(3.5, unit.angstrom),
               unit.Quantity(2.5, unit.angstrom),
               unit.Quantity(2.5, unit.angstrom),
               unit.Quantity(2.5, unit.angstrom),
               unit.Quantity(2.5, unit.angstrom)],
    'epsilons': [unit.Quantity(0.066, unit.kilocalorie / unit.mole),
                 unit.Quantity(0.03, unit.kilocalorie / unit.mole),
                 unit.Quantity(0.03, unit.kilocalorie / unit.mole),
                 unit.Quantity(0.03, unit.kilocalorie / unit.mole),
                 unit.Quantity(0.03, unit.kilocalorie / unit.mole)],
    'bonds': [{'atom1_idx': 0, 'atom2_idx': 1,
               'spring_constant': unit.Quantity(
                   340.0, unit.kilocalorie / (unit.angstrom ** 2
                                              * unit.mole)),
               'eq_dist': unit.Quantity(1.09, unit.angstrom)},
              {'atom1_idx': 0, 'atom2_idx': 2,
               'spring_constant': unit.Quantity(
                   340.0, unit.kilocalorie / (unit.angstrom ** 2
                                              * unit.mole)),
               'eq_dist': unit.Quantity(1.09, unit.angstrom)},
              {'atom1_idx': 0,
               'atom2_idx': 3,
               'spring_constant': unit.Quantity(340.0, unit.kilocalorie
                                                / (unit.angstrom ** 2
                                                   * unit.mole)),
               'eq_dist': unit.Quantity(1.09, unit.angstrom)},
              {'atom1_idx': 0, 'atom2_idx': 4,
               'spring_constant': unit.Quantity(340.0, unit.kilocalorie
                                                / (unit.angstrom ** 2
                                                   * unit.mole)),
               'eq_dist': unit.Quantity(1.09, unit.angstrom)}],
    'angles': [{'atom1_idx': 1, 'atom2_idx': 0, 'atom3_idx': 2,
                'spring_constant': unit.Quantity(
                    33.0, unit.kilocalorie / (unit.mole
                                              * unit.radian ** 2)),
                'eq_angle': unit.Quantity(107.8, unit.degree)},
               {'atom1_idx': 1, 'atom2_idx': 0, 'atom3_idx': 3,
                'spring_constant': unit.Quantity(
                    33.0, unit.kilocalorie / (unit.mole
                                              * unit.radian ** 2)),
                'eq_angle': unit.Quantity(107.8, unit.degree)},
               {'atom1_idx': 1, 'atom2_idx': 0, 'atom3_idx': 4,
                'spring_constant': unit.Quantity(
                    33.0, unit.kilocalorie / (unit.mole
                                              * unit.radian ** 2)),
                'eq_angle': unit.Quantity(107.8, unit.degree)},
               {'atom1_idx': 2, 'atom2_idx': 0, 'atom3_idx': 3,
                'spring_constant': unit.Quantity(
                    33.0, unit.kilocalorie / (unit.mole
                                              * unit.radian ** 2)),
                'eq_angle': unit.Quantity(107.8, unit.degree)},
               {'atom1_idx': 2, 'atom2_idx': 0, 'atom3_idx': 4,
                'spring_constant': unit.Quantity(
                    33.0, unit.kilocalorie / (unit.mole
                                              * unit.radian ** 2)),
                'eq_angle': unit.Quantity(107.8, unit.degree)},
               {'atom1_idx': 3, 'atom2_idx': 0, 'atom3_idx': 4,
                'spring_constant': unit.Quantity(
                    33.0, unit.kilocalorie / (unit.mole
                                              * unit.radian ** 2)),
                'eq_angle': unit.Quantity(107.8, unit.degree)}],
    'SGB_radii': [unit.Quantity(1.975, unit.angstrom),
                  unit.Quantity(1.425, unit.angstrom),
                  unit.Quantity(1.425, unit.angstrom),
                  unit.Quantity(1.425, unit.angstrom),
                  unit.Quantity(1.425, unit.angstrom)],
    'vdW_radii': [unit.Quantity(1.75, unit.angstrom),
                  unit.Quantity(1.25, unit.angstrom),
                  unit.Quantity(1.25, unit.angstrom),
                  unit.Quantity(1.25, unit.angstrom),
                  unit.Quantity(1.25, unit.angstrom)],
    'gammas': [0.005, 0.00859824, 0.00859824, 0.00859824, 0.00859824],
    'alphas': [-0.74168571, 0.268726247, 0.268726247, 0.268726247,
               0.268726247]})


class TestSchrodingerToolkitWrapper(object):
    """
    It wraps all tests that check the SchrodingerToolkitWrapperMolecularGraph
    class.
    """

    def test_get_Schrodinger_parameters(self):
        """
        It tests the standard methods to obtain Schrodinger parameters
        from an offpele's Molecule.
        """

        # Load benzene ring
        molecule = Molecule(smiles='c1ccccc1')
        molecule.parameterize(FORCEFIELD_NAME, charges_method='gasteiger')

        with pytest.raises(ToolkitUnavailableException):
            molecule.get_OPLS_parameters()

        molecule._OPLS_parameters = METHANE_OPLS_PARAMETERS
        _ = molecule.get_OPLS_parameters()

        molecule.add_OPLS_nonbonding_params()
        molecule.add_OPLS_bonds_and_angles()

    def test_assign_Schrodinger_parameters(self):
        """
        It tests the assignment of Schrodinger parameters to offpele's
        Molecule.
        """

        def check_parameters():
            """ It checks the current parameters of the molecule. """
            for atom in molecule.atoms:
                w_atom = WritableAtom(atom)
                w_parameters = [w_atom.index, w_atom.parent.index, w_atom.core,
                                w_atom.OPLS_type, w_atom.PDB_name,
                                w_atom.unknown, w_atom.sigma, w_atom.epsilon,
                                w_atom.charge, w_atom.born_radius,
                                w_atom.SASA_radius, w_atom.nonpolar_gamma,
                                w_atom.nonpolar_alpha]
                assert w_parameters in expected_nonbonding, \
                    'Invalid writable nonbonding parameters {}'.format(w_parameters)

            for bond in molecule.bonds:
                w_bond = WritableBond(bond)
                w_parameters = [attr[1] for attr in list(w_bond)]
                assert w_parameters in expected_bonds, \
                    'Invalid writable bond parameters {}'.format(w_parameters)

            for angle in molecule.angles:
                w_angle = WritableAngle(angle)
                w_parameters = [attr[1] for attr in list(w_angle)]
                assert w_parameters in expected_angles, \
                    'Invalid writable angle parameters {}'.format(w_parameters)

        # Load benzene ring
        molecule = Molecule(smiles='C')
        molecule.parameterize(FORCEFIELD_NAME)

        # Load OPLS parameters
        molecule._OPLS_parameters = METHANE_OPLS_PARAMETERS

        # Check parameters
        # First check
        expected_nonbonding = [
            [1, 0, 'M', 'OFFT', '_C1_', 0, 3.3996695084235347, 0.1094,
             -0.1088, 0, 1.6998347542117673, 0, 0],
            [2, 1, 'M', 'OFFT', '_H1_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [3, 1, 'M', 'OFFT', '_H2_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [4, 1, 'M', 'OFFT', '_H3_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0],
            [5, 1, 'M', 'OFFT', '_H4_', 0, 2.649532787749369, 0.0157,
             0.0267, 0, 1.3247663938746845, 0, 0]]

        expected_bonds = [
            [1, 2, 376.8940758588, 1.094223427522],
            [1, 3, 376.8940758588, 1.094223427522],
            [1, 4, 376.8940758588, 1.094223427522],
            [1, 5, 376.8940758588, 1.094223427522]]

        expected_angles = [
            [2, 1, 3, 33.78875634641, 110.2468561538],
            [2, 1, 4, 33.78875634641, 110.2468561538],
            [2, 1, 5, 33.78875634641, 110.2468561538],
            [3, 1, 4, 33.78875634641, 110.2468561538],
            [3, 1, 5, 33.78875634641, 110.2468561538],
            [4, 1, 5, 33.78875634641, 110.2468561538]]

        check_parameters()

        # Second check
        molecule.add_OPLS_nonbonding_params()

        expected_nonbonding = [
            [1, 0, 'M', 'CT', '_C1_', 0, 3.5, 0.066, -0.1088, 1.975, 1.75,
             0.005, -0.74168571],
            [2, 1, 'M', 'HC', '_H1_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247],
            [3, 1, 'M', 'HC', '_H2_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247],
            [4, 1, 'M', 'HC', '_H3_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247],
            [5, 1, 'M', 'HC', '_H4_', 0, 2.5, 0.03, 0.0267, 1.425, 1.25,
             0.00859824, 0.268726247]]

        check_parameters()

        # Third check
        molecule.add_OPLS_bonds_and_angles()

        expected_bonds = [
            [1, 2, 340.0, 1.09],
            [1, 3, 340.0, 1.09],
            [1, 4, 340.0, 1.09],
            [1, 5, 340.0, 1.09]]

        expected_angles = [
            [2, 1, 3, 33.0, 107.8],
            [2, 1, 4, 33.0, 107.8],
            [2, 1, 5, 33.0, 107.8],
            [3, 1, 4, 33.0, 107.8],
            [3, 1, 5, 33.0, 107.8],
            [4, 1, 5, 33.0, 107.8]]

        check_parameters()

    def test_add_solvent_parameters(self):
        """
        It tests the function that adds the solvent parameters to
        the OPLSParameters collection.
        """

        # Set dummy environment variable
        os.environ["SCHRODINGER"] = "/"

        schrodingerToolkitWrapper = SchrodingerToolkitWrapper()

        # Using a standard atom type
        OPLS_params1 = schrodingerToolkitWrapper.OPLSParameters(
            {'atom_names': [' C1 ', ' H1 ', ' H2 ', ' H3 ', ' H4 '],
             'atom_types': ['CT', 'HC', 'HC', 'HC', 'HC'],
             'charges': [-0.24, 0.06, 0.06, 0.06, 0.06],
             'sigmas': [3.5, 2.5, 2.5, 2.5, 2.5],
             'epsilons': [0.066, 0.03, 0.03, 0.03, 0.03]})

        # Using a similar atom type
        OPLS_params2 = schrodingerToolkitWrapper.OPLSParameters(
            {'atom_names': [' C1 ', ' H1 ', ' H2 ', ' H3 ', ' H4 '],
             'atom_types': ['C3M', 'HC', 'HC', 'HC', 'HC'],
             'charges': [-0.24, 0.06, 0.06, 0.06, 0.06],
             'sigmas': [3.5, 2.5, 2.5, 2.5, 2.5],
             'epsilons': [0.066, 0.03, 0.03, 0.03, 0.03]})

        # Using a default atom type
        OPLS_params3 = schrodingerToolkitWrapper.OPLSParameters(
            {'atom_names': [' C1 ', ' H1 ', ' H2 ', ' H3 ', ' H4 '],
             'atom_types': ['XX', 'HC', 'HC', 'HC', 'HC'],
             'charges': [-0.24, 0.06, 0.06, 0.06, 0.06],
             'sigmas': [3.5, 2.5, 2.5, 2.5, 2.5],
             'epsilons': [0.066, 0.03, 0.03, 0.03, 0.03]})

        schrodingerToolkitWrapper._add_solvent_parameters(OPLS_params1)
        schrodingerToolkitWrapper._add_solvent_parameters(OPLS_params2)
        schrodingerToolkitWrapper._add_solvent_parameters(OPLS_params3)

        assert OPLS_params1['SGB_radii'][0] == \
            unit.Quantity(1.975, unit.angstrom), 'Unexpected SGB radius'
        assert OPLS_params1['vdW_radii'][0] == \
            unit.Quantity(1.750, unit.angstrom), 'Unexpected vdW radius'
        assert OPLS_params1['gammas'][0] == 0.005000000, 'Unexpected gamma'
        assert OPLS_params1['alphas'][0] == -0.741685710, 'Unexpected alpha'

        assert OPLS_params2['SGB_radii'][0] == \
            unit.Quantity(2.002, unit.angstrom), 'Unexpected SGB radius'
        assert OPLS_params2['vdW_radii'][0] == \
            unit.Quantity(1.775, unit.angstrom), 'Unexpected vdW radius'
        assert OPLS_params2['gammas'][0] == 0.023028004, 'Unexpected gamma'
        assert OPLS_params2['alphas'][0] == -0.852763146, 'Unexpected alpha'

        assert OPLS_params3['SGB_radii'][0] == \
            unit.Quantity(1.500, unit.angstrom), 'Unexpected SGB radius'
        assert OPLS_params3['vdW_radii'][0] == \
            unit.Quantity(1.250, unit.angstrom), 'Unexpected vdW radius'
        assert OPLS_params3['gammas'][0] == 0.005000000, 'Unexpected gamma'
        assert OPLS_params3['alphas'][0] == 0.000000000, 'Unexpected alpha'
