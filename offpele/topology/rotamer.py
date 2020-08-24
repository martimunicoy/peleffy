"""
This module handles all classes and functions related with rotamers.
"""

from collections import defaultdict
import networkx as nx
from copy import deepcopy

from offpele.utils.toolkits import RDKitToolkitWrapper


class Rotamer(object):
    """
    It represents the conformational space of a rotatable bond, discretized
    with a certain resolution.
    """

    def __init__(self, index1, index2, resolution=30):
        """
        It initiates an Rotamer object.

        Parameters
        ----------
        atom1 : int
            The index of the first atom involved in the rotamer
        atom2 : int
            The index of the second atom involved in the rotamer
        resolution : float
            The resolution to discretize the rotamer's conformational space
        """
        self._index1 = index1
        self._index2 = index2
        self._resolution = resolution

    @property
    def index1(self):
        """
        Rotamer's atom1 index.

        Returns
        -------
        index1 : int
            The index of the first atom involved in this Rotamer object
        """
        return self._index1

    @property
    def index2(self):
        """
        Rotamer's atom2 index.

        Returns
        -------
        index2 : int
            The index of the second atom involved in this Rotamer object
        """
        return self._index2

    @property
    def resolution(self):
        """
        Rotamer's resolution.

        Returns
        -------
        resolution : float
            The resolution of this Rotamer object
        """
        return self._resolution


class RotamerLibrary(object):
    """
    It represents a set of rotamers found in the same molecule.
    """

    def __init__(self, molecule):
        """
        It initiates a RotamerLibrary object.

        Parameters
        ----------
        molecule : An offpele.topology.Molecule
            The Molecule object whose rotamer library will be generated
        """
        self._molecule = molecule

    def to_file(self, path):
        """
        It writes this RotamerLibrary to a file.

        Parameters
        ----------
        path : str
            Path to save the RotamerLibrary to
        """
        # PELE needs underscores instead of whitespaces
        pdb_atom_names = [name.replace(' ', '_',)
                          for name in self.molecule.get_pdb_atom_names()]

        with open(path, 'w') as file:
            file.write('rot assign res {} &\n'.format(self.molecule.name))
            for i, rotamer_branches in enumerate(self.molecule.rotamers):
                if i > 0:
                    file.write('     newgrp &\n')
                for rotamer in rotamer_branches:
                    atom_name1 = pdb_atom_names[rotamer.index1]
                    atom_name2 = pdb_atom_names[rotamer.index2]
                    file.write('   sidelib FREE{} {} {} &\n'.format(
                        rotamer.resolution, atom_name1, atom_name2))

    @property
    def molecule(self):
        """
        The offpele's Molecule.

        Returns
        -------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        """
        return self._molecule


class MolecularGraph(nx.Graph):
    """
    It represents the structure of a Molecule as a networkx.Graph.
    """

    def __init__(self, molecule):
        """
        It initializes a MolecularGraph object.

        Parameters
        ----------
        molecule : An offpele.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        super().__init__(self)
        self._molecule = molecule
        self._compute_rotamer_graph()
        self._build_core_nodes()

    def _compute_rotamer_graph(self):
        """
        It initializes the netwrokx.Graph with a Molecule object.
        """
        rdkit_toolkit = RDKitToolkitWrapper()
        rot_bonds_atom_ids = \
            rdkit_toolkit.get_atom_ids_with_rotatable_bonds(self.molecule)

        rdkit_molecule = self.molecule.rdkit_molecule

        atom_names = rdkit_toolkit.get_atom_names(self.molecule)

        assert len(atom_names) == len(rdkit_molecule.GetAtoms()), \
            'The length of atom names must match the length of ' \
            + 'molecule\'s atoms'

        for atom, name in zip(rdkit_molecule.GetAtoms(), atom_names):
            self.add_node(atom.GetIdx(), pdb_name=name,
                          nrot_neighbors=list())

        for bond in rdkit_molecule.GetBonds():
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()

            if (frozenset([atom1, atom2]) in rot_bonds_atom_ids):
                rotatable = True
            else:
                rotatable = False
                self.nodes[atom1]['nrot_neighbors'].append(atom2)
                self.nodes[atom2]['nrot_neighbors'].append(atom1)

            self.add_edge(bond.GetBeginAtomIdx(),
                          bond.GetEndAtomIdx(),
                          weight=int(rotatable))

        for i, j in rot_bonds_atom_ids:
            self[i][j]['weight'] = 1
            self.nodes[i]['rotatable'] = True
            self.nodes[j]['rotatable'] = True

    def _build_core_nodes(self):
        """
        It builds the list of core nodes
        """

        def get_all_nrot_neighbors(self, atom_id, visited_neighbors):
            """
            A recursive function that hierarchically visits all atom neighbors
            in the graph.

            Parameters
            ----------
            atom_id : int
                Is is both the id of the graph's node and index of the
                corresponding atom
            visited_neighbors : set[int]
                The ids of the nodes that have already been visited

            Returns
            -------
            visited_neighbors : set[int]
                The updated set that contains the ids of the nodes that have
                already been visited
            """
            if atom_id in visited_neighbors:
                return visited_neighbors
            visited_neighbors.add(atom_id)
            nrot_neighbors = self.nodes[atom_id]['nrot_neighbors']
            for nrot_neighbor in nrot_neighbors:
                visited_neighbors = get_all_nrot_neighbors(
                    self, nrot_neighbor, visited_neighbors)
            return visited_neighbors

        from networkx.algorithms.shortest_paths.generic import \
            shortest_path_length
        from networkx.algorithms.distance_measures import eccentricity

        # Calculate graph distances according to weight values
        weighted_distances = dict(shortest_path_length(self, weight="weight"))

        # Calculate eccentricites using weighted distances
        eccentricities = eccentricity(self, sp=weighted_distances)

        # Group nodes by eccentricity
        nodes_by_eccentricities = defaultdict(list)
        for node, ecc in eccentricities.items():
            nodes_by_eccentricities[ecc].append(node)

        # Core atoms must have the minimum eccentricity
        _, centered_nodes = sorted(nodes_by_eccentricities.items())[0]

        # Construct nrot groups with centered nodes
        # already_visited = set()
        centered_node_groups = list()
        for node in centered_nodes:
            # if node in already_visited:
            #    continue
            centered_node_groups.append(get_all_nrot_neighbors(self, node,
                                                               set()))

        # In case of more than one group, core will be the largest
        core_nodes = sorted(centered_node_groups, key=len, reverse=True)[0]

        # To do: think on what to do with the code below
        """
        # Core can hold a maximum of one rotatable bond <- Not true!
        # Get all core's neighbors
        neighbor_candidates = set()
        for node in core_nodes:
            neighbors = self.neighbors(node)
            for neighbor in neighbors:
                if neighbor not in core_nodes:
                    neighbor_candidates.add(neighbor)

        # If any core's neighbor, get the deepest one and include it to
        # the core
        if len(neighbor_candidates) > 0:
            branch_graph = deepcopy(self)

            for node in core_nodes:
                branch_graph.remove_node(node)

            branch_groups = list(nx.connected_components(branch_graph))

            rot_bonds_per_group = self._get_rot_bonds_per_group(branch_groups)

            best_group = sorted(rot_bonds_per_group, key=len,
                                reverse=True)[0]

            for neighbor in neighbor_candidates:
                if any([neighbor in rot_bond for rot_bond in best_group]):
                    deepest_neighbor = neighbor
                    break
            else:
                raise Exception('Unconsistent graph')

            deepest_neighbors = get_all_nrot_neighbors(self, deepest_neighbor,
                                                       set())

            for neighbor in deepest_neighbors:
                core_nodes.add(neighbor)
        """

        self._core_nodes = core_nodes

    def set_core(self):
        """
        It sets the core of the molecule to minimize the amount of consecutive
        rotamers as much as possible.

        Please, note that the molecule needs to be already parameterized with
        the Open Force Field toolkit before calling this function.
        """
        self.molecule.assert_parameterized()

        for atom in self.molecule.atoms:
            if atom.index in self.core_nodes:
                atom.set_as_core()
            else:
                atom.set_as_branch()

    def set_parents(self):
        """
        It sets the parent of each atom according to the molecular graph.

        Please, note that the molecule needs to be already parameterized with
        the Open Force Field toolkit before calling this function.
        """

        def recursive_child_visitor(parent, already_visited=set()):
            """
            A recursive function that hierarchically visits all the childs of
            each atom.

            Parameters
            ----------
            parent : an offpele.topology.Atom
                The atom whose childs will be visited
            already_visited : set[offpele.topology.Atom]
                The Atom objects that have already been visited

            Returns
            -------
            visited_neighbors : set[offpele.topology.Atom]
                The updated set that contains the Atom objects that have
                already been visited
            """
            if parent in already_visited:
                return already_visited

            already_visited.add(parent)

            childs = [self.molecule.atoms[n] for n in self.neighbors(parent.index)]

            for child in childs:
                if child in already_visited:
                    continue
                child.set_parent(parent)
                already_visited = recursive_child_visitor(child,
                                                          already_visited)

            return already_visited

        self.molecule.assert_parameterized()

        # Start from an atom from the core
        parent = None
        for atom in self.molecule.atoms:
            if atom.core:
                parent = atom
                break

        # Assert a parent was found
        assert parent is not None, 'No core atom found in molecule ' \
            '{}'.format(self.molecule.name)

        already_visited = recursive_child_visitor(parent)

        # Assert all nodes were explored
        assert len(already_visited) == len(self.molecule.atoms), 'Not all ' \
            'nodes were explored'

        # Assert absolut parent is the only with a None parent value
        assert parent.parent is None and \
            sum([int(a.parent is not None) for a in self.molecule.atoms]) \
            == len(self.molecule.atoms) - 1, 'Found descendant without parent'

    def _get_rot_bonds_per_group(self, branch_groups):
        """
        It constructs the rotatable bonds of each branch group.

        Parameters
        ----------
        branch_groups : list[list[int]]
            The node ids of each branch

        Returns
        -------
        rot_bonds_per_group : list[tuple[int, int]]
            The atom ids of all the graph's edges that belong to a rotatable
            bond
        """
        rot_bonds_per_group = list()
        for group in branch_groups:
            rot_bonds = list()
            visited_bonds = set()
            for node in group:
                bonds = self.edges(node)
                for bond in bonds:
                    if bond in visited_bonds:
                        continue
                    if self[bond[0]][bond[1]]['weight'] == 1:
                        rot_bonds.append(bond)
                    visited_bonds.add(bond)
                    visited_bonds.add((bond[1], bond[0]))
            rot_bonds_per_group.append(rot_bonds)

        return rot_bonds_per_group

    def _get_core_atom_per_group(self, rot_bonds_per_group):
        """
        It obtains the core atom for each group.

        Parameters
        ----------
        rot_bonds_per_group : list[tuple[int, int]]
            The atom ids of all the graph's edges that belong to a rotatable
            bond

        Returns
        -------
        core_atom_per_group : list[int]
            The atom id of the atom that belongs to the core for each branch
        """
        core_atom_per_group = list()
        for rot_bonds in rot_bonds_per_group:
            for (a1, a2) in rot_bonds:
                if a1 in self.core_nodes:
                    core_atom_per_group.append(a1)
                    break
                elif a2 in self.core_nodes:
                    core_atom_per_group.append(a2)
                    break
            else:
                core_atom_per_group.append(None)

        return core_atom_per_group

    def _get_sorted_bonds_per_group(self, core_atom_per_group,
                                    rot_bonds_per_group, distances):
        """
        It sorts in increasing order the rotamers of each group according
        to their distance with respect to the corresponding core atom.

        Parameters
        ----------
        core_atom_per_group : list[int]
            The atom id of the atom that belongs to the core for each branch
        rot_bonds_per_group : list[tuple[int, int]]
            The atom ids of all the graph's edges that belong to a rotatable
            bond
        distances : dict[int, dict[int, int]]
            The distance between each pair of nodes (or atoms)

        Returns
        -------
        sorted_rot_bonds_per_group : list[list]
            The rotatable bonds per group, sorted in increasing order by
            their distance with respect to the corresponding core atom
        """
        sorted_rot_bonds_per_group = list()
        for core_atom, rot_bonds in zip(core_atom_per_group,
                                        rot_bonds_per_group):
            sorting_dict = dict()
            for bond in rot_bonds:
                min_d = min([distances[core_atom][bond[0]],
                             distances[core_atom][bond[1]]])
                sorting_dict[bond] = min_d

            sorted_rot_bonds_per_group.append(
                [i[0] for i in
                 sorted(sorting_dict.items(), key=lambda item: item[1])])

        return sorted_rot_bonds_per_group

    def _ignore_terminal_rotatable_bonds(self, sorted_rot_bonds_per_group,
                                         distances):
        """
        It ignores a certain number of terminal rotatable bonds of each
        group.

        Parameters
        ----------
        sorted_rot_bonds_per_group : list[list]
            The rotatable bonds per group, sorted in increasing order by
            their distance with respect to the corresponding core atom
        distances : dict[int, dict[int, int]]
            The distance between each pair of nodes (or atoms)

        Returns
        -------
        filtered_rot_bonds_per_group : list[list]
            The filtered rotatable bonds per group that are obtained
        """
        filtered_rot_bonds_per_group = list()

        # To determine the outter atom in a rotamer, we only need to
        # calculate their distance to any core atom
        core_atom = list(self.core_nodes)[0]

        for rot_bonds in sorted_rot_bonds_per_group:
            rotamer_to_evaluate = rot_bonds[-1]

            node1, node2 = rotamer_to_evaluate

            distance1 = distances[core_atom][node1]
            distance2 = distances[core_atom][node2]

            if distance1 < distance2:
                inner_node = node1
                outter_node = node2
            else:
                inner_node = node2
                outter_node = node1

            # The condition for the current rotamer to be ignored is that the
            # outter node is only attached to terminal nodes (with degree 1)
            ignore_this_rotamer = True
            for neighbor in self.neighbors(outter_node):
                if neighbor is inner_node:
                    continue

                if self.degree(neighbor) > 1:
                    ignore_this_rotamer = False
                    break

            if ignore_this_rotamer:
                filtered_rot_bonds_per_group.append(rot_bonds[:-1])
            else:
                filtered_rot_bonds_per_group.append(rot_bonds)

        return filtered_rot_bonds_per_group

    def get_rotamers(self):
        """
        It builds the RotamerLibrary object.

        Returns
        -------
        rotamers : list[list]
            The list of rotamers grouped by the branch they belong to
        """
        resolution = self.molecule.rotamer_resolution

        assert len(self.core_nodes) > 0, 'No core nodes were found'

        branch_graph = deepcopy(self)

        for node in self.core_nodes:
            branch_graph.remove_node(node)

        branch_groups = list(nx.connected_components(branch_graph))

        rot_bonds_per_group = self._get_rot_bonds_per_group(branch_groups)

        core_atom_per_group = self._get_core_atom_per_group(
            rot_bonds_per_group)

        distances = dict(nx.shortest_path_length(self))

        sorted_rot_bonds_per_group = self._get_sorted_bonds_per_group(
            core_atom_per_group, rot_bonds_per_group, distances)

        """
        if not include_terminal_rotamers:
            sorted_rot_bonds_per_group = \
                self._ignore_terminal_rotatable_bonds(
                    sorted_rot_bonds_per_group, distances)
        """

        rotamers = list()

        # TODO extend core by including one rotatable bond from the largest
        # branch to increase the performance of the algorithm.
        for group_id, rot_bonds in enumerate(sorted_rot_bonds_per_group):
            branch_rotamers = list()
            for (atom1_index, atom2_index) in rot_bonds:
                rotamer = Rotamer(atom1_index, atom2_index, resolution)
                branch_rotamers.append(rotamer)

            if len(branch_rotamers) > 0:
                rotamers.append(branch_rotamers)

        return rotamers

    @property
    def molecule(self):
        """
        The offpele's Molecule.

        Returns
        -------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        """
        return self._molecule

    @property
    def core_nodes(self):
        """
        The list of core nodes.

        Returns
        -------
        core_nodes : list[int]
            The nodes in the core
        """
        return self._core_nodes
