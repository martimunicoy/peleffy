"""
This module contains the tests to check the peleffy's toolkits.
"""

import pytest

from peleffy.topology import Molecule
from peleffy.utils.toolkits import ToolkitUnavailableException
from peleffy.forcefield import OPLS2005ParameterWrapper

from simtk import unit


FORCEFIELD_NAME = 'openff_unconstrained-1.2.0.offxml'
METHANE_OPLS_PARAMETERS = OPLS2005ParameterWrapper({
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
        from an peleffy's Molecule.
        """

        # Load benzene ring
        molecule = Molecule(smiles='c1ccccc1')

        with pytest.raises(ToolkitUnavailableException):
            molecule.parameterize('OPLS2005', charge_method='gasteiger')


class TestRDKitToolkitWrapper(object):
    """
    It wraps all tests that check the RDKitToolkitWrapper class.
    """

    def test_conformer_setter(self):
        """It checks the conformer setter of the RDKit toolkit"""

        from rdkit import Chem
        from copy import deepcopy

        from peleffy.utils import get_data_file_path

        # Load molecule
        mol = Molecule(get_data_file_path('ligands/propionic_acid.pdb'))
        mol.parameterize('openff_unconstrained-1.2.1.offxml',
                         charge_method='gasteiger')

        # Choose a dihedral to track
        dihedral = (0, 1, 2, 3)

        # Get initial dihedral's theta
        conformer = mol.rdkit_molecule.GetConformer()
        initial_theta = Chem.rdMolTransforms.GetDihedralDeg(conformer,
                                                            *dihedral)

        if initial_theta < -179:
            initial_theta += 180.0
        elif initial_theta > 179:
            initial_theta -= 180.0

        assert abs(initial_theta - -0.002) < 10e-3, \
            'Unexpected initial theta value'

        # Get a copy of the rdkit's molecule representation
        rdkit_mol = deepcopy(mol.rdkit_molecule)

        # Modify its conformer
        conformer = rdkit_mol.GetConformer()
        Chem.rdMolTransforms.SetDihedralDeg(conformer, *dihedral, 90)
        new_theta = Chem.rdMolTransforms.GetDihedralDeg(conformer,
                                                        *dihedral)

        assert abs(new_theta - 89.999) < 10e-3, \
            'Unexpected new theta value'

        # Set new conformer to peleffy molecule
        mol.set_conformer(conformer)

        # Check new set theta value
        conformer = mol.rdkit_molecule.GetConformer()
        new_set_theta = Chem.rdMolTransforms.GetDihedralDeg(conformer,
                                                            *dihedral)

        assert abs(new_set_theta - 89.999) < 10e-3, \
            'Unexpected new set theta value'
