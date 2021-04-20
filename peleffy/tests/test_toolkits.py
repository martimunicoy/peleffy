"""
This module contains the tests to check the peleffy's toolkits.
"""

import numpy as np
import pytest

from peleffy.forcefield.parameters import OPLS2005ParameterWrapper
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
        from peleffy.topology import Molecule
        from peleffy.forcefield import OPLS2005ForceField
        from peleffy.utils.toolkits import ToolkitUnavailableException

        # Load benzene ring
        molecule = Molecule(smiles='c1ccccc1')

        # Load OPLS2005 force field
        opls2005 = OPLS2005ForceField()

        # Ensure SCHRODINGER is not in the environment
        import os
        if 'SCHRODINGER' in os.environ:
            del(os.environ['SCHRODINGER'])

        with pytest.raises(ToolkitUnavailableException):
            opls2005.parameterize(molecule, charge_method='gasteiger')


class TestRDKitToolkitWrapper(object):
    """
    It wraps all tests that check the RDKitToolkitWrapper class.
    """

    def test_conformer_setter(self):
        """It checks the conformer setter of the RDKit toolkit"""
        from peleffy.topology import Molecule
        from rdkit import Chem
        from copy import deepcopy

        from peleffy.utils import get_data_file_path

        # Load molecule
        mol = Molecule(get_data_file_path('ligands/propionic_acid.pdb'))

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

    def test_atom_degrees(self):
        """It checks that the atom degree getter works well."""

        from peleffy.topology import Molecule
        from peleffy.utils.toolkits import RDKitToolkitWrapper
        from peleffy.utils import get_data_file_path

        wrapper = RDKitToolkitWrapper()

        pdb_path = get_data_file_path('ligands/methane.pdb')
        m = Molecule(pdb_path)
        degree_by_name = dict(zip(wrapper.get_atom_names(m),
                                  wrapper.get_atom_degrees(m)))

        assert degree_by_name == {' C1 ': 4, ' H1 ': 1, ' H2 ': 1,
                                  ' H3 ': 1, ' H4 ': 1}, \
            'Unexpected pairing between atom names and degrees'

        pdb_path = get_data_file_path('ligands/ethylene.pdb')
        m = Molecule(pdb_path)
        degree_by_name = dict(zip(wrapper.get_atom_names(m),
                                  wrapper.get_atom_degrees(m)))

        assert degree_by_name == {' C1 ': 3, ' C2 ': 3, ' H1 ': 1,
                                  ' H2 ': 1, ' H3 ': 1, ' H4 ': 1}, \
            'Unexpected pairing between atom names and degrees'

        pdb_path = get_data_file_path('ligands/acetylene.pdb')
        m = Molecule(pdb_path)
        degree_by_name = dict(zip(wrapper.get_atom_names(m),
                                  wrapper.get_atom_degrees(m)))

        assert degree_by_name == {' C1 ': 2, ' C2 ': 2, ' H1 ': 1,
                                  ' H2 ': 1}, \
            'Unexpected pairing between atom names and degrees'

        pdb_path = get_data_file_path('ligands/propionic_acid.pdb')
        m = Molecule(pdb_path)
        degree_by_name = dict(zip(wrapper.get_atom_names(m),
                                  wrapper.get_atom_degrees(m)))

        assert degree_by_name == {' C1 ': 4, ' C2 ': 4, ' C3 ': 3,
                                  ' O1 ': 1, ' O2 ': 2, ' H1 ': 1,
                                  ' H2 ': 1, ' H3 ': 1, ' H4 ': 1,
                                  ' H5 ': 1, ' H6 ': 1}, \
            'Unexpected pairing between atom names and degrees'

        pdb_path = get_data_file_path('ligands/trimethylglycine.pdb')
        m = Molecule(pdb_path)
        degree_by_name = dict(zip(wrapper.get_atom_names(m),
                                  wrapper.get_atom_degrees(m)))

        assert degree_by_name == {' C1 ': 4, ' N1 ': 4, ' C2 ': 4,
                                  ' C3 ': 4, ' C4 ': 4, ' C5 ': 3,
                                  ' O1 ': 1, ' O2 ': 1, ' H1 ': 1,
                                  ' H2 ': 1, ' H3 ': 1, ' H4 ': 1,
                                  ' H5 ': 1, ' H6 ': 1, ' H7 ': 1,
                                  ' H8 ': 1, ' H9 ': 1, ' H10': 1,
                                  ' H11': 1}, \
            'Unexpected pairing between atom names and degrees'

        pdb_path = get_data_file_path('ligands/malonate.pdb')
        m = Molecule(pdb_path)
        degree_by_name = dict(zip(wrapper.get_atom_names(m),
                                  wrapper.get_atom_degrees(m)))

        assert degree_by_name == {' O1 ': 1, ' C1 ': 3, ' O2 ': 1,
                                  ' C2 ': 4, ' C3 ': 3, ' O3 ': 2,
                                  ' O4 ': 1, ' H1 ': 1, ' H2 ': 1,
                                  ' H3 ': 1}, \
            'Unexpected pairing between atom names and degrees'

    def test_pdb_parsers(self):
        """It checks that the PDB parsers from RDKit are correctly working."""

        from rdkit.Chem.rdmolfiles import MolToPDBBlock
        from peleffy.utils.toolkits import RDKitToolkitWrapper
        from peleffy.utils import get_data_file_path

        wrapper = RDKitToolkitWrapper()

        pdb_path = get_data_file_path('ligands/benzene.pdb')
        with open(pdb_path) as f:
            pdb_block = f.read()

        rdkit_mol1 = wrapper.from_pdb(pdb_path)
        rdkit_mol2 = wrapper.from_pdb_block(pdb_block)

        block1 = MolToPDBBlock(rdkit_mol1)
        block2 = MolToPDBBlock(rdkit_mol2)

        assert block1 == block2, 'Unexpected pair of RDKit molecules'

    def test_dihedral_angle(self):
        """It checks that the dihedral angle calculator works well."""

        from peleffy.topology import Molecule
        from peleffy.utils.toolkits import RDKitToolkitWrapper
        from peleffy.utils import get_data_file_path

        wrapper = RDKitToolkitWrapper()

        pdb_path = get_data_file_path('ligands/trimethylglycine.pdb')
        m = Molecule(pdb_path)
        dihedral_degrees = wrapper.get_dihedral(m, 2, 1, 4, 5, units="degrees")
        dihedral_rad = wrapper.get_dihedral(m, 2, 1, 4, 5)
        np.testing.assert_almost_equal(dihedral_degrees, -176.348, decimal=2)
        np.testing.assert_almost_equal(dihedral_degrees, np.rad2deg(dihedral_rad), decimal=3)

    def test_dihedral_angle_2(self):
        """It checks that the dihedral angle calculator works well."""

        from peleffy.topology import Molecule
        from peleffy.utils.toolkits import RDKitToolkitWrapper
        from peleffy.utils import get_data_file_path

        wrapper = RDKitToolkitWrapper()

        pdb_path = get_data_file_path('ligands/trimethylglycine.pdb')
        m = Molecule(pdb_path)
        dihedral_degrees = wrapper.get_dihedral(m, 17, 4, 5, 6, units="degrees")
        dihedral_rad = wrapper.get_dihedral(m, 17, 4, 5, 6)
        np.testing.assert_almost_equal(dihedral_degrees, 54.828, decimal=2)
        np.testing.assert_almost_equal(dihedral_degrees, np.rad2deg(dihedral_rad), decimal=3)

    def test_rmsd(self):
        """It checks that the rmsd calculator works well."""

        from peleffy.topology import Molecule
        from peleffy.utils.toolkits import RDKitToolkitWrapper
        from peleffy.utils import get_data_file_path

        wrapper = RDKitToolkitWrapper()

        pdb_path = get_data_file_path('ligands/trimethylglycine.pdb')
        m = Molecule(pdb_path)
        pdb_path2 = get_data_file_path('ligands/trimethylglycine_moved.pdb')
        m2 = Molecule(pdb_path2)
        np.testing.assert_almost_equal(wrapper.get_rmsd(m, m2), 0.3346, decimal=3)
