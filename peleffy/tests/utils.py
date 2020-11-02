"""
This module contains a variety of helpful tools for tests.
"""

import numpy as np
from simtk import unit


SET_OF_LIGAND_PATHS = ['ligands/BNZ.pdb', 'ligands/TOL.pdb', 'ligands/MDB.pdb',
                       'ligands/BIA.pdb', 'ligands/SBN.pdb', 'ligands/OLC.pdb']


def apply_PELE_dihedral_equation(proper, x):
    """
    Given an x, it applies the PELE's dihedral equation to obtain a y.

    Parameters
    ----------
    proper : an peleffy.topology.Proper
        The proper whose parameters will be applied to equation
    x : float
        Equation's x value
    """
    return proper.constant * (1 + proper.prefactor
                              * np.cos(proper.periodicity * x))


def apply_OFF_dihedral_equation(proper, x):
    """
    Given an x, it applies the Open Force Field's dihedral equation to obtain
    a y.

    Parameters
    ----------
    proper : an peleffy.topology.Proper
        The proper whose parameters will be applied to equation
    x : float
        Equation's x value
    """
    return proper.k * (1 + np.cos(proper.periodicity * x - proper.phase))


def check_CHO_charges_in_molecule(molecule):
    """
    It checks whether C, H and O atoms in a molecule have reasonable
    partial charges assigned.

    Parameters
    ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object to check

    Raises
    ------
    AssertionError in case that the checking fails
    ValueError if an unexpected element is found in the molecule
    """
    for atom in molecule.atoms:
        name = atom.PDB_name
        charge = atom.charge.value_in_unit(unit.elementary_charge)

        if 'C' in name:
            assert charge < 1.0 and charge > -0.23, \
                'Presumably not a reasonable charge for a C: ' \
                + '{}'.format(charge)
        elif 'H' in name:
            assert charge > 0 and charge < 0.2, \
                'Presumable not a reasonable charge for a H: ' \
                + '{}'.format(charge)
        elif 'O' in name:
            assert charge < -0.5 and charge > -1.0, \
                'Presumable not a reasonable charge for a O: ' \
                + '{}'.format(charge)
        else:
            raise ValueError('Unknown atom name')


def check_parameters(molecule, expected_nonbonding=None,
                     expected_bonds=None, expected_angles=None,
                     expected_propers=None, expected_impropers=None):
    """
    It checks the current parameters of the molecule.

    Parameters
    ----------
    molecule : an peleffy.topology.Molecule
        The peleffy's Molecule object
    expected_nonbonding : list[list]
        The list of expected nonbonding parameters
    expected_bonds : list[list]
        The list of expected bond parameters
    expected_angles : list[list]
        The list of expected angle parameters
    expected_propers : list[list]
        The list of expected proper parameters
    expected_impropers : list[list]
        The list of expected improper parameters
    """

    from peleffy.template.impact import (WritableAtom,
                                         WritableBond,
                                         WritableAngle,
                                         WritableProper,
                                         WritableImproper)

    if expected_nonbonding is not None:
        assert len(molecule.atoms) == len(expected_nonbonding), \
            'Invalid number of nonbonding terms'

        for atom in molecule.atoms:
            w_atom = WritableAtom(atom)
            w_parameters = [w_atom.index, w_atom.parent.index, w_atom.core,
                            w_atom.OPLS_type, w_atom.PDB_name,
                            w_atom.unknown, w_atom.sigma, w_atom.epsilon,
                            w_atom.charge, w_atom.born_radius,
                            w_atom.SASA_radius, w_atom.nonpolar_gamma,
                            w_atom.nonpolar_alpha]
            assert w_parameters in expected_nonbonding, \
                'Invalid writable nonbonding parameters ' \
                + '{}'.format(w_parameters)

    if expected_bonds is not None:
        assert len(molecule.bonds) == len(expected_bonds), \
            'Invalid number of bond terms'

        for bond in molecule.bonds:
            w_bond = WritableBond(bond)
            w_parameters = [attr[1] for attr in list(w_bond)]
            assert w_parameters in expected_bonds, \
                'Invalid writable bond parameters {}'.format(w_parameters)

    if expected_angles is not None:
        assert len(molecule.angles) == len(expected_angles), \
            'Invalid number of angle terms'

        for angle in molecule.angles:
            w_angle = WritableAngle(angle)
            w_parameters = [attr[1] for attr in list(w_angle)]
            assert w_parameters in expected_angles, \
                'Invalid writable angle parameters {}'.format(w_parameters)

    if expected_propers is not None:
        assert len(molecule.propers) == len(expected_propers), \
            'Invalid number of proper terms'

        for proper in molecule.propers:
            w_proper = WritableProper(proper)
            w_parameters = [attr[1] for attr in list(w_proper)]
            assert w_parameters in expected_propers, \
                'Invalid writable proper parameters {}'.format(w_parameters)

    if expected_impropers is not None:
        assert len(molecule.impropers) == len(expected_impropers), \
            'Invalid number of improper terms'

        for improper in molecule.impropers:
            w_improper = WritableImproper(improper)
            w_parameters = [attr[1] for attr in list(w_improper)]
            assert w_parameters in expected_impropers, \
                'Invalid writable improper parameters {}'.format(w_parameters)


def compare_files(file1, file2):
    """
    It compares two files line by line and gives an AssertionError if there
    is any difference.

    Parameters
    ----------
    file1 : str
        Path to the first file to compare
    file2 : str
        Path to the second file to compare

    Raises
    ------
        AssertionError
            If any difference is found between files
    """
    with open(file1, 'r') as f1:
        lines1 = [line for line in f1.readlines()
                  if not line.startswith("*")]
    with open(file2, 'r') as f2:
        lines2 = [line for line in f2.readlines()
                  if not line.startswith("*")]
    assert len(lines1) == len(lines2), \
        'Number of lines do not match: ' \
        + str(len(lines1)) + ' and ' + str(len(lines2))
    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        assert line1 == line2, \
            'Found different lines at line {}:'.format(i) + '\n' + line1 + line2
