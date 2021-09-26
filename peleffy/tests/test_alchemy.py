"""
This module contains tests that check that the alchemy module.
"""

import pytest


def generate_molecules_and_topologies_from_pdb(pdb1, pdb2):
    """
    It generates the molecules and topologies from two PDB files.

    Parameters
    ----------
    pdb1 : str
        The path to the first PDB
    pdb2 : str
        The path to the second PDB

    Returns
    -------
    molecule1: a peleffy.topology.Molecule
        The first molecule to map
    molecule2: a peleffy.topology.Molecule
        The second molecule to map
    topology1 : a peleffy.topology.Topology object
        The molecular topology representation of molecule 1
    topology2 : a peleffy.topology.Topology object
        The molecular topology representation of molecule 2
    """
    from peleffy.topology import Molecule, Topology
    from peleffy.forcefield import OpenForceField
    from peleffy.utils import get_data_file_path

    mol1 = Molecule(get_data_file_path(pdb1))
    mol2 = Molecule(get_data_file_path(pdb2))

    openff = OpenForceField('openff_unconstrained-2.0.0.offxml')

    params1 = openff.parameterize(mol1, charge_method='gasteiger')
    params2 = openff.parameterize(mol2, charge_method='gasteiger')

    top1 = Topology(mol1, params1)
    top2 = Topology(mol2, params2)

    return mol1, mol2, top1, top2


def generate_molecules_and_topologies_from_smiles(smiles1, smiles2):
    """
    It generates the molecules and topologies from two PDB files.

    Parameters
    ----------
    smiles1 : str
        The SMILES tag of the first molecule
    smiles2 : str
        The SMILES tag of the second molecule

    Returns
    -------
    molecule1: a peleffy.topology.Molecule
        The first molecule to map
    molecule2: a peleffy.topology.Molecule
        The second molecule to map
    topology1 : a peleffy.topology.Topology object
        The molecular topology representation of molecule 1
    topology2 : a peleffy.topology.Topology object
        The molecular topology representation of molecule 2
    """
    from peleffy.topology import Molecule, Topology
    from peleffy.forcefield import OpenForceField

    mol1 = Molecule(smiles=smiles1, hydrogens_are_explicit=False)
    mol2 = Molecule(smiles=smiles2, hydrogens_are_explicit=False)

    openff = OpenForceField('openff_unconstrained-2.0.0.offxml')

    params1 = openff.parameterize(mol1, charge_method='gasteiger')
    params2 = openff.parameterize(mol2, charge_method='gasteiger')

    top1 = Topology(mol1, params1)
    top2 = Topology(mol2, params2)

    return mol1, mol2, top1, top2


class TestAlchemy(object):
    """Alchemy test."""

    def test_alchemizer_initialization_checker(self):
        """
        It checks the initialization checker of Alchemizer class.
        """
        from peleffy.topology import Alchemizer

        mol1, mol2, top1, top2 = \
            generate_molecules_and_topologies_from_smiles('C=C',
                                                          'C(Cl)(Cl)(Cl)')

        _ = Alchemizer(top1, top2)

        with pytest.raises(TypeError):
            _ = Alchemizer(mol1, top2)

        with pytest.raises(TypeError):
            _ = Alchemizer(top1, mol2)

    @pytest.mark.parametrize("pdb1, pdb2, smiles1, smiles2, mapping, " +
                             "non_native_atoms, non_native_bonds, "
                             "non_native_angles, non_native_propers, "
                             "non_native_impropers, exclusive_atoms, "
                             "exclusive_bonds, exclusive_angles, "
                             "exclusive_propers, exclusive_impropers",
                             [(None,
                               None,
                               'C=C',
                               'C(Cl)(Cl)(Cl)',
                               [(0, 0), (1, 1), (2, 2), (3, 4)],
                               [6, ],
                               [],
                               [6, 7, 8],
                               [],
                               [],
                               [4, 5],
                               [3, 4],
                               [0, 1, 5],
                               [],
                               [0, ]
                               ),
                              (None,
                               None,
                               'c1ccccc1',
                               'c1ccccc1C',
                               [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4),
                                (5, 5), (11, 6), (10, 11), (9, 10), (8, 9),
                                (7, 8), (6, 7)],
                               [12, 13, 14],
                               [],
                               [18, 19, 20, 21, 22, 23],
                               [],
                               [],
                               [],
                               [],
                               [],
                               [],
                               []
                               ),
                              (None,
                               None,
                               'c1ccccc1C',
                               'c1ccccc1',
                               [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4),
                                (5, 5), (6, 11), (11, 10), (10, 9), (9, 8),
                                (8, 7), (7, 6)],
                               [],
                               [],
                               [],
                               [],
                               [],
                               [12, 13, 14],
                               [12, 13, 14],
                               [18, 19, 20, 21, 22, 23],
                               [],
                               []
                               ),
                              ('ligands/acetylene.pdb',
                               'ligands/ethylene.pdb',
                               None,
                               None,
                               [(0, 0), (1, 1), (3, 4), (2, 2)],
                               [4, 5],
                               [],
                               [2, 3, 4, 5],
                               [1, 2],
                               [],
                               [],
                               [],
                               [],
                               [],
                               []
                               ),
                              ('ligands/malonate.pdb',
                               'ligands/propionic_acid.pdb',
                               None,
                               None,
                               [(0, 5), (1, 0), (2, 6), (3, 1), (4, 2),
                                (5, 4), (9, 10), (6, 3), (7, 8), (8, 9)],
                               [10],
                               [],
                               [13, 14, 15],
                               [23, 24],
                               [],
                               [],
                               [],
                               [],
                               [],
                               []
                               ),
                              ('ligands/trimethylglycine.pdb',
                               'ligands/propionic_acid.pdb',
                               None,
                               None,
                               [(0, 0), (1, 1), (4, 2), (17, 3), (5, 4),
                                (6, 10), (2, 8), (3, 9), (8, 5), (9, 6),
                                (10, 7)],
                               [],
                               [],
                               [],
                               [],
                               [],
                               [7, 11, 12, 13, 14, 15, 16, 18],
                               [7, 8, 9, 10, 11, 12, 15, 17],
                               [6, 7, 8, 9, 10, 11, 14, 19, 21, 22, 26, 27,
                                28, 29, 30, 31, 32],
                               [40, 41],
                               []
                               ),
                              ('ligands/trimethylglycine.pdb',
                               'ligands/benzamidine.pdb',
                               None,
                               None,
                               [(1, 6), (0, 7), (8, 14), (9, 15), (2, 8),
                                (11, 16)],
                               [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
                               [],
                               [33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,
                                44, 45, 46, 47, 48, 49, 50, 51, 52],
                               [46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
                                57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
                                68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78],
                               [1, 2, 3, 4, 5, 6, 7],
                               [3, 4, 5, 6, 7, 10, 12, 13, 14, 15, 16, 17, 18],
                               [3, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17],
                               [1, 2, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                16, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27,
                                28, 29, 30, 31, 32],
                               [3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17,
                                18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                                29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
                                39, 40, 41, 42, 43, 44, 45],
                               [0]
                               )
                              ])
    def test_alchemizer_initialization(self, pdb1, pdb2, smiles1, smiles2,
                                       mapping, non_native_atoms,
                                       non_native_bonds, non_native_angles,
                                       non_native_propers,
                                       non_native_impropers,
                                       exclusive_atoms, exclusive_bonds,
                                       exclusive_angles, exclusive_propers,
                                       exclusive_impropers):
        """
        It checks the initialization of Alchemizer class.
        """
        from peleffy.topology import Alchemizer

        if pdb1 is not None and pdb2 is not None:
            mol1, mol2, top1, top2 = \
                generate_molecules_and_topologies_from_pdb(pdb1, pdb2)
        elif smiles1 is not None and smiles2 is not None:
            mol1, mol2, top1, top2 = \
                generate_molecules_and_topologies_from_smiles(smiles1, smiles2)
        else:
            raise ValueError('Invalid input parameters for the test')

        alchemizer = Alchemizer(top1, top2)

        # Check alchemizer content
        assert alchemizer._mapping == mapping, \
            'Unexpected mapping'
        assert alchemizer._non_native_atoms == non_native_atoms, \
            'Unexpected non native atoms'
        assert alchemizer._non_native_bonds == non_native_bonds, \
            'Unexpected non native bonds'
        assert alchemizer._non_native_angles == non_native_angles, \
            'Unexpected non native angles'
        assert alchemizer._non_native_propers == non_native_propers, \
            'Unexpected non native propers'
        assert alchemizer._non_native_impropers == non_native_impropers, \
            'Unexpected non native impropers'
        assert alchemizer._exclusive_atoms == exclusive_atoms, \
            'Unexpected exclusive atoms'
        assert alchemizer._exclusive_bonds == exclusive_bonds, \
            'Unexpected exclusive bonds'
        assert alchemizer._exclusive_angles == exclusive_angles, \
            'Unexpected exclusive angles'
        assert alchemizer._exclusive_propers == exclusive_propers, \
            'Unexpected exclusive propers'
        assert alchemizer._exclusive_impropers == exclusive_impropers, \
            'Unexpected exclusive impropers'

    def test_fep_lambda(self):
        """
        It validates the effects of fep lambda on atom parameters.
        """
        from peleffy.topology import Alchemizer
        from peleffy.template.impact import (WritableAtom, WritableBond,
                                             WritableAngle, WritableProper,
                                             WritableImproper)

        mol1, mol2, top1, top2 = \
            generate_molecules_and_topologies_from_smiles('C=C',
                                                          'C(Cl)(Cl)(Cl)')

        alchemizer = Alchemizer(top1, top2)

        top = alchemizer.get_alchemical_topology(fep_lambda=0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert (sigma2 / sigma1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert (epsilon2 / epsilon1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert (SASA_radius2 / SASA_radius1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert (charge2 / charge1) - (1- 0.2) < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert (bond_sc2 / bond_sc1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert (angle_sc2 / angle_sc1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert (proper_c2 / proper_c1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert (improper_c2 / improper_c1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=1.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.4)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert (sigma2 / sigma1) - 0.4 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert (epsilon2 / epsilon1) - 0.4 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert (SASA_radius2 / SASA_radius1) - 0.4 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert (charge2 / charge1) - 0.4 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert (bond_sc2 / bond_sc1) - 0.4 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert (angle_sc2 / angle_sc1) - 0.4 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert (proper_c2 / proper_c1) - 0.4 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert (improper_c2 / improper_c1) - 0.4 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas1.append(atom.sigma)
                epsilons1.append(atom.epsilon)
                SASA_radii1.append(atom.SASA_radius)
                charges1.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants1.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants1.append(improper.constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=1.0)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas2.append(atom.sigma)
                epsilons2.append(atom.epsilon)
                SASA_radii2.append(atom.SASA_radius)
                charges2.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants2.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants2.append(improper.constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

    def test_coul1_lambda(self):
        """
        It validates the effects of coul1 lambda on atom parameters.
        """
        from peleffy.topology import Alchemizer
        from peleffy.template.impact import (WritableAtom, WritableBond,
                                             WritableAngle, WritableProper,
                                             WritableImproper)

        mol1, mol2, top1, top2 = \
            generate_molecules_and_topologies_from_smiles('C=C',
                                                          'C(Cl)(Cl)(Cl)')

        alchemizer = Alchemizer(top1, top2)

        top = alchemizer.get_alchemical_topology(fep_lambda=0,
                                                 coul1_lambda=0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul1_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert (charge2 / charge1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul1_lambda=0.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul1_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul1_lambda=0.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas1.append(atom.sigma)
                epsilons1.append(atom.epsilon)
                SASA_radii1.append(atom.SASA_radius)
                charges1.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants1.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants1.append(improper.constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=1.0,
                                                 coul1_lambda=1.0)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas2.append(atom.sigma)
                epsilons2.append(atom.epsilon)
                SASA_radii2.append(atom.SASA_radius)
                charges2.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants2.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants2.append(improper.constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

    def test_coul2_lambda(self):
        """
        It validates the effects of coul2 lambda on atom parameters.
        """
        from peleffy.topology import Alchemizer
        from peleffy.template.impact import (WritableAtom, WritableBond,
                                             WritableAngle, WritableProper,
                                             WritableImproper)

        mol1, mol2, top1, top2 = \
            generate_molecules_and_topologies_from_smiles('C=C',
                                                          'C(Cl)(Cl)(Cl)')

        alchemizer = Alchemizer(top1, top2)

        top = alchemizer.get_alchemical_topology(fep_lambda=0,
                                                 coul2_lambda=0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul2_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul2_lambda=1.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul2_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert (charge2 / charge1) - 0.8 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 coul2_lambda=0.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas1.append(atom.sigma)
                epsilons1.append(atom.epsilon)
                SASA_radii1.append(atom.SASA_radius)
                charges1.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants1.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants1.append(improper.constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=1.0,
                                                 coul2_lambda=1.0)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas2.append(atom.sigma)
                epsilons2.append(atom.epsilon)
                SASA_radii2.append(atom.SASA_radius)
                charges2.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants2.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants2.append(improper.constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

    def test_vdw_lambda(self):
        """
        It validates the effects of vdw lambda on atom parameters.
        """
        from peleffy.topology import Alchemizer
        from peleffy.template.impact import (WritableAtom, WritableBond,
                                             WritableAngle, WritableProper,
                                             WritableImproper)

        mol1, mol2, top1, top2 = \
            generate_molecules_and_topologies_from_smiles('C=C',
                                                          'C(Cl)(Cl)(Cl)')

        alchemizer = Alchemizer(top1, top2)

        top = alchemizer.get_alchemical_topology(fep_lambda=0,
                                                 vdw_lambda=0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 vdw_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._exclusive_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._exclusive_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._exclusive_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._exclusive_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._exclusive_propers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert (sigma2 / sigma1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert (epsilon2 / epsilon1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert (SASA_radius2 / SASA_radius1) - (1 - 0.2) < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 vdw_lambda=1.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas1.append(atom.sigma)
            epsilons1.append(atom.epsilon)
            SASA_radii1.append(atom.SASA_radius)
            charges1.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants1.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants1.append(improper.spring_constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 vdw_lambda=0.2)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in alchemizer._non_native_atoms:
            atom = WritableAtom(top.atoms[atom_idx])
            sigmas2.append(atom.sigma)
            epsilons2.append(atom.epsilon)
            SASA_radii2.append(atom.SASA_radius)
            charges2.append(atom.charge)

        for bond_idx in alchemizer._non_native_bonds:
            bond = WritableBond(top.bonds[bond_idx])
            bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in alchemizer._non_native_angles:
            angle = WritableAngle(top.angles[angle_idx])
            angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in alchemizer._non_native_propers:
            proper = WritableProper(top.propers[proper_idx])
            proper_constants2.append(proper.spring_constant)

        for improper_idx in alchemizer._non_native_impropers:
            improper = WritableImproper(top.impropers[improper_idx])
            improper_constants2.append(improper.spring_constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert (sigma2 / sigma1) - 0.2 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert (epsilon2 / epsilon1) - 0.2 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert (SASA_radius2 / SASA_radius1) - 0.2 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

        top = alchemizer.get_alchemical_topology(fep_lambda=0.0,
                                                 vdw_lambda=0.0)

        sigmas1 = list()
        epsilons1 = list()
        SASA_radii1 = list()
        charges1 = list()
        bond_spring_constants1 = list()
        angle_spring_constants1 = list()
        proper_constants1 = list()
        improper_constants1 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas1.append(atom.sigma)
                epsilons1.append(atom.epsilon)
                SASA_radii1.append(atom.SASA_radius)
                charges1.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants1.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants1.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants1.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants1.append(improper.constant)

        top = alchemizer.get_alchemical_topology(fep_lambda=1.0,
                                                 vdw_lambda=1.0)

        sigmas2 = list()
        epsilons2 = list()
        SASA_radii2 = list()
        charges2 = list()
        bond_spring_constants2 = list()
        angle_spring_constants2 = list()
        proper_constants2 = list()
        improper_constants2 = list()

        for atom_idx in range(0, len(top.atoms)):
            if (atom_idx not in alchemizer._exclusive_atoms and
                    atom_idx not in alchemizer._non_native_atoms):
                atom = WritableAtom(top.atoms[atom_idx])
                sigmas2.append(atom.sigma)
                epsilons2.append(atom.epsilon)
                SASA_radii2.append(atom.SASA_radius)
                charges2.append(atom.charge)

        for bond_idx in range(0, len(top.bonds)):
            if (bond_idx not in alchemizer._exclusive_bonds and
                    bond_idx not in alchemizer._non_native_bonds):
                bond = WritableBond(top.bonds[bond_idx])
                bond_spring_constants2.append(bond.spring_constant)

        for angle_idx in range(0, len(top.angles)):
            if (angle_idx not in alchemizer._exclusive_angles and
                    angle_idx not in alchemizer._non_native_angles):
                angle = WritableAngle(top.angles[angle_idx])
                angle_spring_constants2.append(angle.spring_constant)

        for proper_idx in range(0, len(top.propers)):
            if (proper_idx not in alchemizer._exclusive_propers and
                    proper_idx not in alchemizer._non_native_propers):
                proper = WritableProper(top.propers[proper_idx])
                proper_constants2.append(proper.constant)

        for improper_idx in range(0, len(top.impropers)):
            if (improper_idx not in alchemizer._exclusive_impropers and
                    improper_idx not in alchemizer._non_native_impropers):
                improper = WritableImproper(top.impropers[improper_idx])
                improper_constants2.append(improper.constant)

        for sigma1, sigma2 in zip(sigmas1, sigmas2):
            assert sigma2 - sigma1 < 1e-5, \
                'Unexpected ratio between sigmas'

        for epsilon1, epsilon2 in zip(epsilons1, epsilons2):
            assert epsilon2 - epsilon1 < 1e-5, \
                'Unexpected ratio between epsilons'

        for SASA_radius1, SASA_radius2 in zip(SASA_radii1, SASA_radii2):
            assert SASA_radius2 - SASA_radius1 < 1e-5, \
                'Unexpected ratio between SASA radii'

        for charge1, charge2 in zip(charges1, charges2):
            assert charge2 - charge1 < 1e-5, \
                'Unexpected ratio between charges'

        for bond_sc1, bond_sc2 in zip(bond_spring_constants1,
                                      bond_spring_constants2):
            assert bond_sc2 - bond_sc1 < 1e-5, \
                'Unexpected ratio between bond spring constants'

        for angle_sc1, angle_sc2 in zip(angle_spring_constants1,
                                        angle_spring_constants2):
            assert angle_sc2 - angle_sc1 < 1e-5, \
                'Unexpected ratio between angle spring constants'

        for proper_c1, proper_c2 in zip(proper_constants1,
                                        proper_constants2):
            assert proper_c2 - proper_c1 < 1e-5, \
                'Unexpected ratio between proper constants'

        for improper_c1, improper_c2 in zip(improper_constants1,
                                            improper_constants2):
            assert improper_c2 - improper_c1 < 1e-5, \
                'Unexpected ratio between improper constants'

    @pytest.mark.parametrize("pdb1, pdb2, smiles1, smiles2, " +
                             "fep_lambda, coul1_lambda, coul2_lambda, " +
                             "vdw_lambda, bonded_lambda, " +
                             "golden_sigmas, golden_epsilons, " +
                             "golden_born_radii, golden_SASA_radii, " +
                             "golden_nonpolar_gammas, " +
                             "golden_nonpolar_alphas",
                             [(None,
                               None,
                               'C=C',
                               'C(Cl)(Cl)(Cl)',
                               0.0,
                               None,
                               None,
                               None,
                               None,
                               [3.480646886945065, 3.480646886945065,
                                2.5725815350632795, 2.5725815350632795,
                                2.5725815350632795, 2.5725815350632795, 0.0],
                               [0.0868793154488, 0.0868793154488,
                                0.01561134320353, 0.01561134320353,
                                0.01561134320353, 0.01561134320353, 0.0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [1.7403234434725325, 1.7403234434725325,
                                1.2862907675316397, 1.2862907675316397,
                                1.2862907675316397, 1.2862907675316397, 0.0],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]
                               ),
                              (None,
                               None,
                               'C=C',
                               'C(Cl)(Cl)(Cl)',
                               0.2,
                               None,
                               None,
                               None,
                               None,
                               [3.480646886945065, 3.480646886945065,
                                2.5725815350632795, 2.5725815350632795,
                                2.0580652280506238, 2.0580652280506238,
                                0.6615055612921249],
                               [0.0868793154488, 0.0868793154488,
                                0.01561134320353, 0.01561134320353,
                                0.012489074562824,
                                0.012489074562824, 0.05312002093054],
                               [0, 0, 0, 0, 0, 0, 0],
                               [1.7403234434725325, 1.7403234434725325,
                                1.2862907675316397, 1.2862907675316397,
                                1.0290326140253119, 1.0290326140253119,
                                0.33075278064606245],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]
                               ),
                              (None,
                               None,
                               'C=C',
                               'C(Cl)(Cl)(Cl)',
                               0.8,
                               None,
                               None,
                               None,
                               None,
                               [3.480646886945065, 3.480646886945065,
                                2.5725815350632795, 2.5725815350632795,
                                0.5145163070126558, 0.5145163070126558,
                                2.6460222451684996],
                               [0.0868793154488, 0.0868793154488,
                                0.01561134320353, 0.01561134320353,
                                0.003122268640705999, 0.003122268640705999,
                                0.21248008372216],
                               [0, 0, 0, 0, 0, 0, 0],
                               [1.7403234434725325, 1.7403234434725325,
                                1.2862907675316397, 1.2862907675316397,
                                0.2572581535063279, 0.2572581535063279,
                                1.3230111225842498],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]
                               ),
                              (None,
                               None,
                               'C=C',
                               'C(Cl)(Cl)(Cl)',
                               1.0,
                               None,
                               None,
                               None,
                               None,
                               [3.480646886945065, 3.480646886945065,
                                2.5725815350632795, 2.5725815350632795,
                                0.0, 0.0,
                                3.3075278064606244],
                               [0.0868793154488, 0.0868793154488,
                                0.01561134320353, 0.01561134320353, 0.0, 0.0,
                                0.2656001046527],
                               [0, 0, 0, 0, 0, 0, 0],
                               [1.7403234434725325, 1.7403234434725325,
                                1.2862907675316397, 1.2862907675316397, 0.0,
                                0.0, 1.6537639032303122],
                               [0, 0, 0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0, 0, 0]
                               ),
                              ])
    def test_atoms_in_alchemical_topology(self, pdb1, pdb2, smiles1, smiles2,
                                          fep_lambda, coul1_lambda,
                                          coul2_lambda, vdw_lambda,
                                          bonded_lambda,
                                          golden_sigmas,
                                          golden_epsilons,
                                          golden_born_radii,
                                          golden_SASA_radii,
                                          golden_nonpolar_gammas,
                                          golden_nonpolar_alphas):
        """
        It validates the effects of lambda on atom parameters.
        """
        from peleffy.topology import Alchemizer
        from peleffy.template.impact import WritableAtom

        mol1, mol2, top1, top2 = \
            generate_molecules_and_topologies_from_smiles(smiles1, smiles2)

        alchemizer = Alchemizer(top1, top2)

        top = alchemizer.get_alchemical_topology(fep_lambda=fep_lambda,
                                                 coul1_lambda=coul1_lambda,
                                                 coul2_lambda=coul2_lambda,
                                                 vdw_lambda=vdw_lambda,
                                                 bonded_lambda=bonded_lambda)

        sigmas = list()
        epsilons = list()
        born_radii = list()
        SASA_radii = list()
        nonpolar_gammas = list()
        nonpolar_alphas = list()

        for atom in top.atoms:
            atom = WritableAtom(atom)
            sigmas.append(atom.sigma)
            epsilons.append(atom.epsilon)
            born_radii.append(atom.born_radius)
            SASA_radii.append(atom.SASA_radius)
            nonpolar_gammas.append(atom.nonpolar_gamma)
            nonpolar_alphas.append(atom.nonpolar_alpha)

        assert sigmas == golden_sigmas, 'Unexpected sigmas'
        assert epsilons == golden_epsilons, 'Unexpected epsilons'
        assert born_radii == golden_born_radii, 'Unexpected born radii'
        assert SASA_radii == golden_SASA_radii, 'Unexpected SASA radii'
        assert nonpolar_gammas == golden_nonpolar_gammas, \
            'Unexpected non polar gammas'
        assert nonpolar_alphas == golden_nonpolar_alphas, \
            'Unexpected non polar alphas'
