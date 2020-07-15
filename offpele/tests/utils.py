"""
This module contains a variety of helpful tools for tests.
"""

import numpy as np


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
