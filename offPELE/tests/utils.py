import numpy as np


SET_OF_LIGAND_PATHS = ['ligands/BNZ.pdb', 'ligands/TOL.pdb', 'ligands/MDB.pdb',
                       'ligands/BIA.pdb', 'ligands/SBN.pdb', 'ligands/OLC.pdb']


def apply_PELE_dihedral_equation(proper, x):
    return proper.constant * (1 + proper.prefactor
                              * np.cos(proper.periodicity * x))


def apply_OFF_dihedral_equation(proper, x):
    return proper.k * (1 + np.cos(proper.periodicity * x - proper.phase))
