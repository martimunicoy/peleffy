"""
This module contains the tests to check the offpele's rotamer library
builder.
"""

import pytest

from offpele.utils import get_data_file_path
from offpele.topology import Molecule
from offpele.topology.rotamer import MolecularGraph


class TestMolecularGraph(object):
    """
    It wraps all tests that check the MolecularGraph class.
    """

    def test_rotamer_library_builder(self):
        """
        It tests the rotamer library builder.
        """
        FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'
        LIGAND_PATH = 'ligands/OLC.pdb'

        ligand_path = get_data_file_path(LIGAND_PATH)
        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME, charges_method='gasteiger')

        graph = MolecularGraph(molecule)

        graph.set_core()
        graph.set_parents()

        rotamer_library = graph.build_rotamer_library(n_rot_bonds_to_ignore=0)

        rotamers_per_branch = rotamer_library.rotamers

        assert len(rotamers_per_branch) == 2, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list_1 = list()
        atom_list_2 = list()
        for branch, rotamers in rotamers_per_branch.items():
            if branch == 0:
                for rotamer in rotamers:
                    atom_list_1.append(set([rotamer.atom1, rotamer.atom2]))
            elif branch == 1:
                for rotamer in rotamers:
                    atom_list_2.append(set([rotamer.atom1, rotamer.atom2]))

        EXPECTED_ATOMS_1 = [set(['_C8_', '_C9_']), set(['_C7_', '_C8_']),
                            set(['_C6_', '_C7_']), set(['_C5_', '_C6_']),
                            set(['_C4_', '_C5_']), set(['_C3_', '_C4_']),
                            set(['_C3_', '_C2_']), set(['_C2_', '_C1_'])]

        EXPECTED_ATOMS_2 = [set(['_C11', '_C10']), set(['_C11', '_C12']),
                            set(['_C12', '_C13']), set(['_C13', '_C14']),
                            set(['_C14', '_C15']), set(['_C15', '_C16']),
                            set(['_C16', '_C17']), set(['_C17', '_C18'])]

        where_1 = list()
        for atom_pair in atom_list_1:
            if atom_pair in EXPECTED_ATOMS_1:
                where_1.append(1)
            elif atom_pair in EXPECTED_ATOMS_2:
                where_1.append(2)
            else:
                where_1.append(0)

        where_2 = list()
        for atom_pair in atom_list_2:
            if atom_pair in EXPECTED_ATOMS_1:
                where_2.append(1)
            elif atom_pair in EXPECTED_ATOMS_2:
                where_2.append(2)
            else:
                where_2.append(0)

        assert (all(i == 1 for i in where_1)
                and all(i == 2 for i in where_2)) or \
            (all(i == 2 for i in where_1)
             and all(i == 1 for i in where_2)), "Invalid rotamer library " + \
            "{}, {}".format(where_1, where_2)

    def test_terminal_rotamer_filtering(self):
        """
        It tests the rotamer library builder when the terminal rotatable bonds
        are ignored.
        """
        FORCEFIELD_NAME = 'openff_unconstrained-1.1.1.offxml'
        LIGAND_PATH = 'ligands/OLC.pdb'

        ligand_path = get_data_file_path(LIGAND_PATH)
        molecule = Molecule(ligand_path)
        molecule.parameterize(FORCEFIELD_NAME, charges_method='gasteiger')

        graph = MolecularGraph(molecule)

        graph.set_core()
        graph.set_parents()

        rotamer_library = graph.build_rotamer_library(n_rot_bonds_to_ignore=1)

        rotamers_per_branch = rotamer_library.rotamers

        assert len(rotamers_per_branch) == 2, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list_1 = list()
        atom_list_2 = list()
        for branch, rotamers in rotamers_per_branch.items():
            if branch == 0:
                for rotamer in rotamers:
                    atom_list_1.append(set([rotamer.atom1, rotamer.atom2]))
            elif branch == 1:
                for rotamer in rotamers:
                    atom_list_2.append(set([rotamer.atom1, rotamer.atom2]))

        EXPECTED_ATOMS_1 = [set(['_C8_', '_C9_']), set(['_C7_', '_C8_']),
                            set(['_C6_', '_C7_']), set(['_C5_', '_C6_']),
                            set(['_C4_', '_C5_']), set(['_C3_', '_C4_']),
                            set(['_C3_', '_C2_'])]

        EXPECTED_ATOMS_2 = [set(['_C11', '_C10']), set(['_C11', '_C12']),
                            set(['_C12', '_C13']), set(['_C13', '_C14']),
                            set(['_C14', '_C15']), set(['_C15', '_C16']),
                            set(['_C16', '_C17'])]

        where_1 = list()
        for atom_pair in atom_list_1:
            if atom_pair in EXPECTED_ATOMS_1:
                where_1.append(1)
            elif atom_pair in EXPECTED_ATOMS_2:
                where_1.append(2)
            else:
                where_1.append(0)

        where_2 = list()
        for atom_pair in atom_list_2:
            if atom_pair in EXPECTED_ATOMS_1:
                where_2.append(1)
            elif atom_pair in EXPECTED_ATOMS_2:
                where_2.append(2)
            else:
                where_2.append(0)

        assert (all(i == 1 for i in where_1)
                and all(i == 2 for i in where_2)) or \
            (all(i == 2 for i in where_1)
             and all(i == 1 for i in where_2)), "Invalid rotamer library " + \
            "{}, {}".format(where_1, where_2)
