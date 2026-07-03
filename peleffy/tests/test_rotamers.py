"""
This module contains the tests to check the peleffy's rotamer library
builder.
"""

import pytest
import networkx as nx
from copy import deepcopy

from peleffy.utils import get_data_file_path
from peleffy.topology import Molecule
from peleffy.topology.rotamer import MolecularGraph


class TestMolecularGraph(object):
    """
    It wraps all tests that check the MolecularGraph class.
    """

    def test_rotamer_library_builder(self):
        """
        It tests the rotamer library builder.
        """
        LIGAND_PATH = 'ligands/oleic_acid.pdb'

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
        LIGAND_PATH = 'ligands/oleic_acid.pdb'

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

        LIGAND_PATH = 'ligands/oleic_acid.pdb'
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

        LIGAND_PATH = 'ligands/oleic_acid.pdb'
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


class TestBifurcationAndRingDetection(object):
    """
    Tests for the ring-detection and branch-splitting logic inside
    MolecularGraph:  _identify_rigid_rings, _split_branch_at_rings,
    and the end-to-end rotamer groups produced for molecules that
    contain rings in their branches.
    """

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _branch_groups(graph):
        """Return the connected components after removing all core nodes."""
        bg = deepcopy(graph)
        for node in graph.core_nodes:
            bg.remove_node(node)
        return list(nx.connected_components(bg))

    @staticmethod
    def _rotamer_sets(molecule):
        """Return each branch as a list of frozensets {idx1, idx2}."""
        return [
            [frozenset([r.index1, r.index2]) for r in branch]
            for branch in molecule.rotamers
        ]

    # ------------------------------------------------------------------
    # _identify_rigid_rings
    # ------------------------------------------------------------------

    def test_identify_rings_linear_molecule(self):
        """A purely linear molecule has no rigid rings in any branch."""
        molecule = Molecule(smiles='CCCCCC', name='hexane', tag='HEX')
        graph = MolecularGraph(molecule)
        for bg in self._branch_groups(graph):
            rings = graph._identify_rigid_rings(bg)
            assert rings == [], \
                "Expected no rigid rings in hexane branch, got {}".format(rings)

    def test_identify_rings_benzene_core(self):
        """
        propylbenzene: the benzene ring is (part of) the core, so the
        branch that contains the ring atoms should still be detected.
        The branch containing the ring atoms (3-8) should yield exactly
        one rigid ring with 6 members.
        """
        molecule = Molecule(smiles='CCCc1ccccc1', name='propylbenzene',
                            tag='PRB')
        graph = MolecularGraph(molecule)
        # atom 2 is the sole core node; branch groups are [0,1,2] (the
        # propyl tail) and [3..8] (the phenyl ring).  Only the ring branch
        # matters here.
        for bg in self._branch_groups(graph):
            rings = graph._identify_rigid_rings(bg)
            if len(bg) == 6:
                # This is the phenyl branch
                assert len(rings) == 1, \
                    "Expected 1 rigid ring in phenyl branch, got {}".format(rings)
                assert len(list(rings)[0]) == 6, \
                    "Expected 6-membered ring, got {}".format(rings)

    def test_identify_rings_branch_contains_ring(self):
        """
        CCCCCc1ccc(CCC)c(CC)c1 — the pentyl chain is one branch; the
        other branch contains the benzene ring together with a propyl and
        an ethyl chain.  _identify_rigid_rings must return exactly one
        6-membered ring for that branch.
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        graph = MolecularGraph(molecule)
        # The ring-containing branch has 11 nodes (indices 5-15)
        ring_branch = None
        for bg in self._branch_groups(graph):
            if len(bg) > 6:
                ring_branch = bg
                break
        assert ring_branch is not None, "Could not find ring-containing branch"
        rings = graph._identify_rigid_rings(ring_branch)
        assert len(rings) == 1, \
            "Expected 1 rigid ring in branch, got {}".format(rings)
        assert len(list(rings)[0]) == 6, \
            "Expected 6-membered ring, got size {}".format(len(list(rings)[0]))

    def test_identify_rings_no_false_positives_for_chains(self):
        """
        In the short-chain branches of a trisubstituted benzene only
        straight chains are present; _identify_rigid_rings must return
        empty lists for those branches.
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        graph = MolecularGraph(molecule)
        for bg in self._branch_groups(graph):
            if len(bg) <= 4:    # the short propyl / ethyl chains
                rings = graph._identify_rigid_rings(bg)
                assert rings == [], \
                    "Unexpected ring in chain branch {}: {}".format(bg, rings)

    # ------------------------------------------------------------------
    # _split_branch_at_rings
    # ------------------------------------------------------------------

    def test_split_no_ring(self):
        """
        A branch with no ring should come back as a single group
        containing all its nodes.
        """
        molecule = Molecule(smiles='CCCCCC', name='hexane', tag='HEX')
        graph = MolecularGraph(molecule)
        for bg in self._branch_groups(graph):
            groups = graph._split_branch_at_rings(bg)
            # All nodes must be covered by exactly one group
            covered = set().union(*groups)
            assert covered == bg, \
                "Split groups do not cover all branch nodes"
            assert len(groups) == 1, \
                "Linear branch should not be split into multiple groups"

    def test_split_non_bifurcating_ring(self):
        """
        propylbenzene — the phenyl ring in the branch has only ONE
        rotatable connection (toward the core), so it is *not* a
        bifurcation point and the branch must remain a single group.
        """
        molecule = Molecule(smiles='CCCc1ccccc1', name='propylbenzene',
                            tag='PRB')
        graph = MolecularGraph(molecule)
        for bg in self._branch_groups(graph):
            groups = graph._split_branch_at_rings(bg)
            assert len(groups) == 1, \
                "Non-bifurcating ring branch should not be split, " \
                "but got {} groups".format(len(groups))

    def test_split_bifurcating_ring_produces_two_groups(self):
        """
        CCCCCc1ccc(CCC)c(CC)c1 — the ring in the large branch has three
        rotatable connections: the pentyl chain (dominant) and two
        minority chains (propyl and ethyl).  _split_branch_at_rings must
        return exactly two groups for that branch.
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        graph = MolecularGraph(molecule)
        ring_branch = max(self._branch_groups(graph), key=len)
        groups = graph._split_branch_at_rings(ring_branch)
        assert len(groups) == 2, \
            "Bifurcating ring should split branch into 2 groups, " \
            "got {}".format(len(groups))

    def test_split_dominant_group_is_larger(self):
        """
        The dominant group (longest chain through the ring) must contain
        more rotatable bonds than any minority group.
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        graph = MolecularGraph(molecule)
        ring_branch = max(self._branch_groups(graph), key=len)
        groups = graph._split_branch_at_rings(ring_branch)

        def rot_bonds_in(nodes):
            count = 0
            seen = set()
            for n in nodes:
                for nb in graph.neighbors(n):
                    if nb in nodes and graph[n][nb]['weight'] == 1:
                        key = frozenset([n, nb])
                        if key not in seen:
                            seen.add(key)
                            count += 1
            return count

        rot_counts = [rot_bonds_in(g) for g in groups]
        assert max(rot_counts) == rot_counts[0], \
            "First group (dominant) should have the most rotatable bonds, " \
            "got counts {}".format(rot_counts)

    # ------------------------------------------------------------------
    # End-to-end rotamer branch tests
    # ------------------------------------------------------------------

    def test_linear_molecule_two_branches(self):
        """
        hexane (CCCCCC): two branches from the central core atom (index 2).
        Branch 0 should contain bonds (2,3) and (3,4); branch 1 bond (2,1).
        """
        molecule = Molecule(smiles='CCCCCC', name='hexane', tag='HEX',
                            exclude_terminal_rotamers=False)
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "hexane should have 2 branches, got {}".format(len(rotamers))
        all_bond_sets = [
            frozenset([r.index1, r.index2])
            for branch in rotamers for r in branch
        ]
        expected = [
            frozenset([2, 3]), frozenset([3, 4]),
            frozenset([2, 1]),
        ]
        for bond in expected:
            assert bond in all_bond_sets, \
                "Expected bond {} missing from hexane rotamers".format(bond)
        assert len(all_bond_sets) == 3, \
            "hexane should have 3 rotatable bonds total, got {}".format(
                len(all_bond_sets))

    def test_ring_core_two_side_chains(self):
        """
        diethylbenzene (CCc1ccc(CC)cc1): benzene ring forms the core;
        two independent ethyl branches.
        """
        molecule = Molecule(smiles='CCc1ccc(CC)cc1', name='diethylbenzene',
                            tag='DEB')
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "diethylbenzene should have 2 branches, got {}".format(
                len(rotamers))
        # Each branch must have exactly 1 rotatable bond
        for i, branch in enumerate(rotamers):
            assert len(branch) == 1, \
                "Each ethyl branch should have 1 rotamer, " \
                "branch {} has {}".format(i, len(branch))

    def test_ring_in_branch_not_bifurcating(self):
        """
        propylbenzene (CCCc1ccccc1): only 2 rotatable bonds (C-C chain);
        the ring branch and the chain branch are separate with no
        bifurcation.  Total branches == 2, total rotatable bonds == 2.
        """
        molecule = Molecule(smiles='CCCc1ccccc1', name='propylbenzene',
                            tag='PRB')
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "propylbenzene should have 2 branches, got {}".format(
                len(rotamers))
        total_rot = sum(len(b) for b in rotamers)
        assert total_rot == 2, \
            "propylbenzene should have 2 rotatable bonds, got {}".format(
                total_rot)

    def test_bifurcating_ring_in_branch_three_groups(self):
        """
        CCCCCc1ccc(CCC)c(CC)c1 (pentyl_trisubst): the benzene ring sits
        inside a branch and has 3 rotatable connections.  After splitting
        there must be 3 rotamer groups: the dominant pentyl-chain group,
        the propyl minority group, and the short ethyl minority group.
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        rotamers = molecule.rotamers
        assert len(rotamers) == 3, \
            "pentyl_trisubst should have 3 rotamer branches, got {}".format(
                len(rotamers))

    def test_bifurcating_ring_dominant_branch_is_longest(self):
        """
        In CCCCCc1ccc(CCC)c(CC)c1 the dominant branch (through the
        pentyl chain) must be the longest group.
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        rotamers = molecule.rotamers
        lengths = [len(b) for b in rotamers]
        assert lengths[0] >= max(lengths), \
            "Dominant branch should be first and longest, got {}".format(lengths)

    def test_bifurcating_ring_bond_membership(self):
        """
        CCCCCc1ccc(CCC)c(CC)c1: verify specific rotatable bonds appear
        in the expected groups.
        Atom ordering (heavy atoms only, 0-indexed):
          0-4  : C-C-C-C-C  (pentyl chain, core at atom 4)
          5    : ring entry from core
          6-15 : phenyl ring + propyl and ethyl substituents
        Expected rotamer groups (sorted by decreasing length):
          - branch 0 (pentyl tail, 3 bonds): (4,3), (3,2), (2,1)
          - branch 1 (dominant ring-side, 3 bonds): (4,5), (8,9), (9,10)
          - branch 2 (minority ethyl, 1 bond): (12,13)
        """
        molecule = Molecule(smiles='CCCCCc1ccc(CCC)c(CC)c1',
                            name='pentyl_trisubst', tag='PTS')
        all_sets = self._rotamer_sets(molecule)

        # branch 1 (ring-side dominant) must include the ring-entry bond (4,5)
        # and the propyl chain bonds (8,9), (9,10)
        dominant_bonds = [frozenset([4, 5]), frozenset([8, 9]), frozenset([9, 10])]
        dominant_flat = [b for b in all_sets[1]]
        for bond in dominant_bonds:
            assert bond in dominant_flat, \
                "Bond {} missing from dominant ring-side branch".format(bond)

        # the minority ethyl bond (12,13) must appear in branch 2
        minority_bond = frozenset([12, 13])
        minority_flat = [b for b in all_sets[2]]
        assert minority_bond in minority_flat, \
            "Bond {} missing from minority ethyl branch".format(minority_bond)

        # branch 0 is the pure pentyl tail with bonds (4,3), (3,2), (2,1)
        tail_bonds = [frozenset([4, 3]), frozenset([3, 2]), frozenset([2, 1])]
        tail_flat = [b for b in all_sets[0]]
        for bond in tail_bonds:
            assert bond in tail_flat, \
                "Bond {} missing from pentyl-tail branch".format(bond)

    def test_bhp_two_symmetric_branches(self):
        """
        BHP.pdb: a 4-hydroxyphenyl molecule with two chains.  The core
        contains the ring; each branch has 3 rotatable bonds.
        """
        ligand_path = get_data_file_path('ligands/BHP.pdb')
        molecule = Molecule(ligand_path)
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "BHP should have 2 branches, got {}".format(len(rotamers))
        # Both branches should have the same number of rotatable bonds (3)
        lengths = sorted([len(b) for b in rotamers])
        assert lengths == [3, 3], \
            "BHP branches should each have 3 rotamers, got {}".format(lengths)
        # Verify known bond pairs (from actual output)
        all_flat = {frozenset([r.index1, r.index2])
                    for branch in rotamers for r in branch}
        expected_bonds = {
            frozenset([6, 0]), frozenset([0, 9]), frozenset([9, 10]),
            frozenset([6, 4]), frozenset([4, 3]), frozenset([3, 5]),
        }
        assert all_flat == expected_bonds, \
            "BHP rotamer bonds mismatch.\n  got:      {}\n  expected: {}".format(
                all_flat, expected_bonds)

    def test_sam_two_branches_with_ring(self):
        """
        SAM.pdb: two branches; one is a short linear chain and the other
        passes through a ring system.  Must produce exactly 2 branches
        with the expected bond pairs.
        """
        ligand_path = get_data_file_path('ligands/SAM.pdb')
        molecule = Molecule(ligand_path)
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "SAM should have 2 branches, got {}".format(len(rotamers))
        all_flat = {frozenset([r.index1, r.index2])
                    for branch in rotamers for r in branch}
        expected_bonds = {
            frozenset([7, 6]), frozenset([6, 5]), frozenset([5, 1]),
            frozenset([1, 2]),
            frozenset([7, 9]), frozenset([9, 10]), frozenset([16, 17]),
        }
        assert all_flat == expected_bonds, \
            "SAM rotamer bonds mismatch.\n  got:      {}\n  expected: {}".format(
                all_flat, expected_bonds)

    def test_lig1_two_branches_ring_in_branch(self):
        """
        LIG1.pdb: two branches each containing 4 rotatable bonds,
        with at least one ring structure inside a branch.
        """
        ligand_path = get_data_file_path('ligands/LIG1.pdb')
        molecule = Molecule(ligand_path)
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "LIG1 should have 2 branches, got {}".format(len(rotamers))
        lengths = sorted([len(b) for b in rotamers])
        assert lengths == [4, 4], \
            "LIG1 branches should each have 4 rotamers, got {}".format(lengths)
        all_flat = {frozenset([r.index1, r.index2])
                    for branch in rotamers for r in branch}
        expected_bonds = {
            frozenset([9, 8]), frozenset([8, 5]), frozenset([5, 4]),
            frozenset([4, 0]),
            frozenset([9, 12]), frozenset([12, 13]), frozenset([13, 14]),
            frozenset([20, 21]),
        }
        assert all_flat == expected_bonds, \
            "LIG1 rotamer bonds mismatch.\n  got:      {}\n  expected: {}".format(
                all_flat, expected_bonds)

    def test_lig2_two_branches(self):
        """
        LIG2.pdb: two branches with 2 and 1 rotatable bonds respectively.
        """
        ligand_path = get_data_file_path('ligands/LIG2.pdb')
        molecule = Molecule(ligand_path)
        rotamers = molecule.rotamers
        assert len(rotamers) == 2, \
            "LIG2 should have 2 branches, got {}".format(len(rotamers))
        all_flat = {frozenset([r.index1, r.index2])
                    for branch in rotamers for r in branch}
        expected_bonds = {
            frozenset([13, 7]), frozenset([7, 6]),
            frozenset([14, 15]),
        }
        assert all_flat == expected_bonds, \
            "LIG2 rotamer bonds mismatch.\n  got:      {}\n  expected: {}".format(
                all_flat, expected_bonds)
