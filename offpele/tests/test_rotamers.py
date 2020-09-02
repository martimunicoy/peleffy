"""
This module contains the tests to check the offpele's rotamer library
builder.
"""

import pytest

from offpele.utils import get_data_file_path
from offpele.topology import Molecule


class TestMolecularGraph(object):
    """
    It wraps all tests that check the MolecularGraph class.
    """

    def test_rotamer_library_builder(self):
        """
        It tests the rotamer library builder.
        """
        LIGAND_PATH = 'ligands/OLC.pdb'

        ligand_path = get_data_file_path(LIGAND_PATH)
        molecule = Molecule(ligand_path, exclude_terminal_rotamers=False)

        # rotamer_library = RotamerLibrary(molecule)

        rotamers_per_branch = molecule.rotamers

        assert len(rotamers_per_branch) == 2, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list_1 = list()
        atom_list_2 = list()

        rotamers = rotamers_per_branch[0]
        for rotamer in rotamers:
            atom_list_1.append(set([rotamer.index1, rotamer.index2]))

        rotamers = rotamers_per_branch[1]
        for rotamer in rotamers:
            atom_list_2.append(set([rotamer.index1, rotamer.index2]))

        EXPECTED_INDICES_1 = [set([9, 10]), set([8, 9]), set([7, 8]),
                              set([6, 7]), set([5, 6]), set([2, 5]),
                              set([0, 2]), set([0, 1])]

        EXPECTED_INDICES_2 = [set([12, 11]), set([12, 13]), set([13, 14]),
                              set([14, 15]), set([15, 16]), set([16, 17]),
                              set([17, 18]), set([18, 19])]

        where_1 = list()
        for atom_pair in atom_list_1:
            if atom_pair in EXPECTED_INDICES_1:
                where_1.append(1)
            elif atom_pair in EXPECTED_INDICES_2:
                where_1.append(2)
            else:
                where_1.append(0)

        where_2 = list()
        for atom_pair in atom_list_2:
            if atom_pair in EXPECTED_INDICES_1:
                where_2.append(1)
            elif atom_pair in EXPECTED_INDICES_2:
                where_2.append(2)
            else:
                where_2.append(0)

        assert (all(i == 1 for i in where_1)
                and all(i == 2 for i in where_2)) or \
            (all(i == 2 for i in where_1)
             and all(i == 1 for i in where_2)), "Invalid rotamer library " + \
            "{}, {}".format(where_1, where_2)

        assert (all(i == 1 for i in where_1)
                and all(i == 2 for i in where_2)
                and len(where_1) == len(EXPECTED_INDICES_1)
                and len(where_2) == len(EXPECTED_INDICES_2)) or \
               (all(i == 2 for i in where_1)
                and all(i == 1 for i in where_2)
                and len(where_1) == len(EXPECTED_INDICES_2)
                and len(where_2) == len(EXPECTED_INDICES_1)), "Unexpected " + \
            "number of rotamers"

    def test_terminal_rotamer_filtering(self):
        """
        It tests the rotamer library builder when the terminal rotatable bonds
        are ignored.
        """
        LIGAND_PATH = 'ligands/OLC.pdb'

        ligand_path = get_data_file_path(LIGAND_PATH)
        molecule = Molecule(ligand_path, exclude_terminal_rotamers=True)

        rotamers_per_branch = molecule.rotamers

        assert len(rotamers_per_branch) == 2, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list_1 = list()
        atom_list_2 = list()
        rotamers = rotamers_per_branch[0]
        for rotamer in rotamers:
            atom_list_1.append(set([rotamer.index1, rotamer.index2]))

        rotamers = rotamers_per_branch[1]
        for rotamer in rotamers:
            atom_list_2.append(set([rotamer.index1, rotamer.index2]))

        EXPECTED_INDICES_1 = [set([9, 10]), set([8, 9]), set([7, 8]),
                              set([6, 7]), set([5, 6]), set([2, 5]),
                              set([0, 2]), set([0, 1])]

        EXPECTED_INDICES_2 = [set([12, 11]), set([12, 13]), set([13, 14]),
                              set([14, 15]), set([15, 16]), set([16, 17]),
                              set([17, 18])]

        where_1 = list()
        for atom_pair in atom_list_1:
            if atom_pair in EXPECTED_INDICES_1:
                where_1.append(1)
            elif atom_pair in EXPECTED_INDICES_2:
                where_1.append(2)
            else:
                where_1.append(0)

        where_2 = list()
        for atom_pair in atom_list_2:
            if atom_pair in EXPECTED_INDICES_1:
                where_2.append(1)
            elif atom_pair in EXPECTED_INDICES_2:
                where_2.append(2)
            else:
                where_2.append(0)

        assert (all(i == 1 for i in where_1)
                and all(i == 2 for i in where_2)) or \
            (all(i == 2 for i in where_1)
             and all(i == 1 for i in where_2)), "Invalid rotamer library " + \
            "{}, {}".format(where_1, where_2)

        assert (all(i == 1 for i in where_1)
                and all(i == 2 for i in where_2)
                and len(where_1) == len(EXPECTED_INDICES_1)
                and len(where_2) == len(EXPECTED_INDICES_2)) or \
               (all(i == 2 for i in where_1)
                and all(i == 1 for i in where_2)
                and len(where_1) == len(EXPECTED_INDICES_2)
                and len(where_2) == len(EXPECTED_INDICES_1)), "Unexpected " + \
            "number of rotamers"
