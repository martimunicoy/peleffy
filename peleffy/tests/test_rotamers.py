"""
This module contains the tests to check the peleffy's rotamer library
builder.
"""

import pytest

from peleffy.utils import get_data_file_path
from peleffy.topology import Molecule


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

    def test_rotamer_core_constraint(self):
        """
        It tests the rotamer library builder when constraining its core
        to contain a specific atom.
        """

        LIGAND_PATH = 'ligands/OLC.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        # Test atom index constraint
        molecule = Molecule(ligand_path, core_constraints=[19, ],
                            exclude_terminal_rotamers=False)

        rotamers_per_branch = molecule.rotamers

        assert len(rotamers_per_branch) == 1, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list = list()
        for rotamer in rotamers_per_branch[0]:
            atom_list.append(set([rotamer.index1, rotamer.index2]))

        EXPECTED_INDICES = [set([18, 19]), set([17, 18]), set([16, 17]),
                            set([15, 16]), set([14, 15]), set([13, 14]),
                            set([12, 13]), set([11, 12]), set([9, 10]),
                            set([8, 9]), set([7, 8]), set([6, 7]),
                            set([5, 6]), set([2, 5]), set([0, 2]),
                            set([0, 1])]

        assert len(atom_list) == len(EXPECTED_INDICES), "Unexpected " + \
            "number of rotamers"

        assert all(atom_pair in EXPECTED_INDICES for atom_pair in atom_list), \
            "Invalid rotamer library"

        # Test PDB atom name constraint
        molecule = Molecule(ligand_path, core_constraints=[' C18', ],
                            exclude_terminal_rotamers=False)

        rotamers_per_branch = molecule.rotamers

        assert len(rotamers_per_branch) == 1, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list = list()
        for rotamer in rotamers_per_branch[0]:
            atom_list.append(set([rotamer.index1, rotamer.index2]))

        EXPECTED_INDICES = [set([18, 19]), set([17, 18]), set([16, 17]),
                            set([15, 16]), set([14, 15]), set([13, 14]),
                            set([12, 13]), set([11, 12]), set([9, 10]),
                            set([8, 9]), set([7, 8]), set([6, 7]),
                            set([5, 6]), set([2, 5]), set([0, 2]),
                            set([0, 1])]

        assert len(atom_list) == len(EXPECTED_INDICES), "Unexpected " + \
            "number of rotamers"

        assert all(atom_pair in EXPECTED_INDICES for atom_pair in atom_list), \
            "Invalid rotamer library"

        # Test core constraint with terminal exclusion
        molecule = Molecule(ligand_path, core_constraints=[' C18', ],
                            exclude_terminal_rotamers=True)

        rotamers_per_branch = molecule.rotamers

        assert len(rotamers_per_branch) == 1, "Found an invalid number " + \
            "of branches: {}".format(len(rotamers_per_branch))

        atom_list = list()
        for rotamer in rotamers_per_branch[0]:
            atom_list.append(set([rotamer.index1, rotamer.index2]))

        EXPECTED_INDICES = [set([17, 18]), set([16, 17]), set([15, 16]),
                            set([14, 15]), set([13, 14]), set([12, 13]),
                            set([11, 12]), set([9, 10]), set([8, 9]),
                            set([7, 8]), set([6, 7]), set([5, 6]),
                            set([2, 5]), set([0, 2]), set([0, 1])]

        assert len(atom_list) == len(EXPECTED_INDICES), "Unexpected " + \
            "number of rotamers"

        assert all(atom_pair in EXPECTED_INDICES for atom_pair in atom_list), \
            "Invalid rotamer library"

        # Test core constraint with a central core
        molecule = Molecule(ligand_path, core_constraints=[' C9 ', ],
                            exclude_terminal_rotamers=True)

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

        # Test core constraint with a multiple central core
        molecule = Molecule(ligand_path,
                            core_constraints=[' C8 ', ' C9 ', ' C10'],
                            exclude_terminal_rotamers=True)

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

        EXPECTED_INDICES_1 = [set([8, 9]), set([7, 8]), set([6, 7]),
                              set([5, 6]), set([2, 5]), set([0, 2]),
                              set([0, 1])]

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

    def test_rotamer_core_constraint_adjacency(self):
        """
        It tests the adjacency check up that is performed prior building
        the rotamer library builder with core constraints.
        """

        LIGAND_PATH = 'ligands/OLC.pdb'
        ligand_path = get_data_file_path(LIGAND_PATH)

        # Test adjacent core constraint selection
        _ = Molecule(ligand_path,
                     core_constraints=[' C8 ', ' C9 ', ' C10'])

        # Test non adjacent core constraint selection
        with pytest.raises(ValueError) as e:
            _ = Molecule(ligand_path,
                         core_constraints=[' C1 ', ' C9 ', ' C10'])

        assert str(e.value) == 'All atoms in atom constraints must be ' \
            + 'adjacent and atom C1 is not'
