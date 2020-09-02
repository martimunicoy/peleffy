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
    proper : an offpele.topology.Proper
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
    proper : an offpele.topology.Proper
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
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object to check

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
