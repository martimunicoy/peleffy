"""
This module contains the tests to check the offpele's toolkits.
"""

import pytest

import os

from offpele.utils.toolkits import SchrodingerToolkitWrapper

from simtk import unit


class TestSchrodingerToolkitWrapper(object):
    """
    It wraps all tests that check the SchrodingerToolkitWrapperMolecularGraph
    class.
    """

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
