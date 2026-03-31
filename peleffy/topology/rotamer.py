"""
This module handles all classes and functions related with rotamers.
"""

from collections import defaultdict
import networkx as nx
from copy import deepcopy

from peleffy.utils.toolkits import RDKitToolkitWrapper


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

    def __eq__(self, other):
        """Define equality operator"""
        return (self.index1 == other.index1
                and self.index2 == other.index2) \
            or (self.index1 == other.index2
                and self.index2 == other.index1)

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
        molecule : An peleffy.topology.Molecule
            The Molecule object whose rotamer library will be generated

        Load a molecule and create its rotamer library template

        >>> from peleffy.topology import Molecule
        >>> from peleffy.topology import RotamerLibrary

        >>> molecule = Molecule(smiles='CCCC', name='butane', tag='BUT')

        >>> rotamer_library = RotamerLibrary(mol)
        >>> rotamer_library.to_file('butz')

        Load a molecule and create its rotamer library template with
        a core constraint

        >>> from peleffy.topology import Molecule
        >>> from peleffy.topology import RotamerLibrary

        >>> molecule = Molecule(smiles='CCCC', name='butane', tag='BUT',
                                exclude_terminal_rotamers=False,
                                core_constraints=[0, ])

        >>> rotamer_library = RotamerLibrary(mol)
        >>> rotamer_library.to_file('butz')

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
            file.write('rot assign res {} &\n'.format(self.molecule.tag))
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
        The peleffy's Molecule.

        Returns
        -------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        return self._molecule

    def _ipython_display_(self):
        """
        It displays a 2D molecular representation with bonds highlighted
        according to this rotamer library object.

        Returns
        -------
        representation_2D : a IPython display object
            It is displayable RDKit molecule with an embeded 2D
            representation
        """
        COLORS = [(82 / 255, 215 / 255, 255 / 255), (255 / 255, 154 / 255, 71 / 255),
                  (161 / 255, 255 / 255, 102 / 255), (255 / 255, 173 / 255, 209 / 255),
                  (154 / 255, 92 / 255, 255 / 255), (66 / 255, 255 / 255, 167 / 255),
                  (251 / 255, 255 / 255, 17 / 255)]

        from IPython.display import display

        # Get 2D molecular representation
        rdkit_toolkit = RDKitToolkitWrapper()
        representation = rdkit_toolkit.get_2D_representation(self.molecule)

        # Get rotamer branches from molecule
        rotamer_branches = self.molecule.rotamers

        bond_indexes = list()
        bond_color_dict = dict()
        for bond in representation.GetBonds():
            rotamer = Rotamer(bond.GetBeginAtom().GetIdx(),
                              bond.GetEndAtom().GetIdx())

            for color_index, group in enumerate(rotamer_branches):
                if rotamer in group:
                    bond_indexes.append(bond.GetIdx())
                    try:
                        bond_color_dict[bond.GetIdx()] = COLORS[color_index]
                    except IndexError:
                        bond_color_dict[bond.GetIdx()] = (99 / 255,
                                                          122 / 255,
                                                          126 / 255)
                    break

        atom_indexes = list()
        radii_dict = dict()
        atom_color_dict = dict()

        for atom in representation.GetAtoms():
            atom_index = atom.GetIdx()
            if atom_index in self.molecule._graph.core_nodes:
                atom_indexes.append(atom.GetIdx())
                radii_dict[atom.GetIdx()] = 0.6
                atom_color_dict[atom.GetIdx()] = (255 / 255, 243 / 255, 133 / 255)

        # Get 2D molecular representation
        rdkit_toolkit = RDKitToolkitWrapper()
        representation = rdkit_toolkit.get_2D_representation(self.molecule)

        # Get its image
        image = rdkit_toolkit.draw_molecule(representation,
                                            atom_indexes=atom_indexes,
                                            radii_dict=radii_dict,
                                            atom_color_dict=atom_color_dict,
                                            bond_indexes=bond_indexes,
                                            bond_color_dict=bond_color_dict)

        return display(image)


class MolecularGraph(nx.Graph):
    """
    It represents the structure of a Molecule as a networkx.Graph.
    """

    def __init__(self, molecule):
        """
        It initializes a MolecularGraph object.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        super().__init__(self)
        self._molecule = molecule
        self._compute_rotamer_graph()
        self._build_core_nodes()

    def _compute_rotamer_graph(self):
        """
        It initializes the network.Graph with a Molecule object.
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

    def get_parents(self, parent):
        """
        It sets the parent of each atom according to the molecular graph.

        Parameters
        ----------
        parent : int
            The index of the node to use as the absolute parent

        Returns
        -------
        parents : dict[int, int]
            A dictionary containing the index of the parent of each
            atom according to the molecular graph, keyed by the index
            of each child
        """

        def recursive_child_visitor(parent, parents,
                                    already_visited=set()):
            """
            A recursive function that hierarchically visits all the childs of
            each atom.

            Parameters
            ----------
            parent : int
                The index of the atom whose childs will be visited
            parents : dict[int, int]
                A dictionary containing the index of the parent of each
                atom according to the molecular graph, keyed by the index
            visited_neighbors : set[int]
                The updated set that contains the indexes of the atoms
                that have already been visited

            Returns
            -------
            parents : dict[int, int]
                A dictionary containing the index of the parent of each
                atom according to the molecular graph, keyed by the index
                of each child
            visited_neighbors : set[int]
                The updated set that contains the indexes of the atoms
                that have already been visited
            """
            if parent in already_visited:
                return already_visited

            already_visited.add(parent)

            childs = self.neighbors(parent)

            for child in childs:
                if child in already_visited:
                    continue
                parents[child] = parent
                parents, already_visited = recursive_child_visitor(
                    child, parents, already_visited)

            return parents, already_visited

        # Initialize the parents dictionary
        parents = {parent: None}

        parents, already_visited = recursive_child_visitor(parent, parents)

        # Assert absolut parent is the only with a None parent value
        if parents[parent] is not None or \
                sum([int(parents[i] is not None) for i in self.nodes]) \
                != len(self.nodes) - 1:

            from peleffy.utils import Logger
            logger = Logger()
            logger.error('Error: found descendant without parent')

        return parents

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
        # Collect all rotatable bonds in the molecule (deduplicated)
        all_rotatable_bonds = []
        seen = set()
        for node in self.nodes():
            for neighbor in self.neighbors(node):
                if self[node][neighbor]['weight'] == 1:
                    key = frozenset([node, neighbor])
                    if key not in seen:
                        seen.add(key)
                        all_rotatable_bonds.append((node, neighbor))

        # Pre-compute distances once (used for tiebreaking below)
        distances = dict(nx.shortest_path_length(self))

        def dist_to_core(atom):
            """
            Minimum graph-distance from atom to any core node.
            
            Parameters
            ----------
            atom : int
                The index of the atom to measure distance from
            
            Returns
            -------
            d : float
                The minimum graph-distance from atom to any core node, or
                infinity if no path exists
            """
            d = float('inf')
            for core_node in self.core_nodes:
                if atom in distances and core_node in distances[atom]:
                    d = min(d, distances[atom][core_node])
            return d

        # For each bond, decide which group it belongs to.
        # Rules (in priority order):
        #   1. If one atom is in the core → assign to the group containing the other atom.
        #   2. If both atoms appear in a common group → assign to the largest such group
        #      (handles shared ring nodes between dominant and minority groups).
        #   3. If atoms belong to different groups → assign to the group whose atom
        #      is closer to core (i.e. the more "core-side" group).
        #   4. Fallback: whichever group contains either atom (largest wins).

        rot_bonds_per_group = [[] for _ in branch_groups]
        assigned_bonds = set()   # frozenset keys of already-assigned bonds

        for bond in all_rotatable_bonds:
            key = frozenset(bond)
            if key in assigned_bonds:
                continue

            atom1, atom2 = bond

            # Which groups contain each atom?
            groups_with_atom1 = [i for i, g in enumerate(branch_groups) if atom1 in g]
            groups_with_atom2 = [i for i, g in enumerate(branch_groups) if atom2 in g]

            in_core1 = atom1 in self.core_nodes
            in_core2 = atom2 in self.core_nodes

            target_group_idx = None

            if in_core1 and groups_with_atom2:
                # atom1 is in the core → the bond belongs to the branch group
                # that contains atom2; if atom2 spans multiple groups, keep
                # the largest one to stay with the dominant chain.
                target_group_idx = max(groups_with_atom2,
                                       key=lambda i: len(branch_groups[i]))

            elif in_core2 and groups_with_atom1:
                # Symmetric case: atom2 is in the core.
                target_group_idx = max(groups_with_atom1,
                                       key=lambda i: len(branch_groups[i]))

            elif groups_with_atom1 and groups_with_atom2:
                common = set(groups_with_atom1) & set(groups_with_atom2)
                if common:
                    # Both atoms appear in the same group(s) – e.g. a ring
                    # node shared between the dominant and a minority group.
                    # Pick the largest group to keep the bond with the dominant
                    # chain.
                    target_group_idx = max(common,
                                           key=lambda i: len(branch_groups[i]))
                else:
                    # Atoms belong to different groups (e.g. a bond that
                    # straddles a ring boundary). Assign the bond to the group
                    # whose atom sits closer to the core so that the kinematic
                    # chain remains directed from core outward.
                    d1 = dist_to_core(atom1)
                    d2 = dist_to_core(atom2)
                    if d1 <= d2:
                        target_group_idx = max(groups_with_atom1,
                                               key=lambda i: len(branch_groups[i]))
                    else:
                        target_group_idx = max(groups_with_atom2,
                                               key=lambda i: len(branch_groups[i]))

            elif groups_with_atom1:
                # Only atom1 appears in a known group; assign the bond there.
                target_group_idx = max(groups_with_atom1,
                                       key=lambda i: len(branch_groups[i]))

            elif groups_with_atom2:
                # Only atom2 appears in a known group; assign the bond there.
                target_group_idx = max(groups_with_atom2,
                                       key=lambda i: len(branch_groups[i]))

            else:
                # Neither atom maps to any branch group – this can happen for
                # bonds that connect exclusively to core nodes. Skip silently.
                continue

            rot_bonds_per_group[target_group_idx].append(bond)
            assigned_bonds.add(key)

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
            core_atom = None
            
            if not rot_bonds:
                core_atom_per_group.append(None)
                continue
            
            # Check if any bond atom is directly in the core
            for (a1, a2) in rot_bonds:
                if a1 in self.core_nodes:
                    core_atom = a1
                    break
                elif a2 in self.core_nodes:
                    core_atom = a2
                    break
            
            # If no direct core atom found, find the closest core atom
            if core_atom is None:
                # Get all atoms in this branch
                branch_atoms = set()
                for bond in rot_bonds:
                    branch_atoms.add(bond[0])
                    branch_atoms.add(bond[1])
                
                # Find the atom in this branch that is closest to any core atom
                min_distance = float('inf')
                distances = dict(nx.shortest_path_length(self))
                
                for atom in branch_atoms:
                    for core_node in self.core_nodes:
                        if atom in distances and core_node in distances[atom]:
                            if distances[atom][core_node] < min_distance:
                                min_distance = distances[atom][core_node]
                                core_atom = core_node
            
            core_atom_per_group.append(core_atom)

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
            if core_atom is None or not rot_bonds:
                sorted_rot_bonds_per_group.append([])
                continue
                
            sorting_dict = dict()
            for bond in rot_bonds:
                if (core_atom in distances and 
                    bond[0] in distances[core_atom] and 
                    bond[1] in distances[core_atom]):
                    if (distances[core_atom][bond[0]] < 
                            distances[core_atom][bond[1]]):
                        sorting_dict[(bond[0], bond[1])] = \
                            distances[core_atom][bond[0]]
                    else:
                        sorting_dict[(bond[1], bond[0])] = \
                            distances[core_atom][bond[1]]
                else:
                    # Fallback: keep original bond order
                    sorting_dict[bond] = 0

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

    def _identify_rigid_rings(self, branch_nodes):
        """
        Identify rigid rings (sub-cores) within a branch using cycle detection.
        
        Parameters
        ----------
        branch_nodes : set[int]
            The nodes in the branch to analyze
            
        Returns
        -------
        rigid_rings : list[set[int]]
            List of node sets, each representing a rigid ring
        """
        # Build a subgraph from all bonds within the branch (rotatable and
        # non-rotatable alike) so that ring-detection works on the full
        # topology rather than just the rigid scaffold.
        subgraph = nx.Graph()
        subgraph.add_nodes_from(branch_nodes)
        
        non_rotatable_bonds = []
        for node in branch_nodes:
            for neighbor in self.neighbors(node):
                if neighbor in branch_nodes:
                    subgraph.add_edge(node, neighbor)
                    if self[node][neighbor]['weight'] == 0:
                        non_rotatable_bonds.append((node, neighbor))
        
        # Enumerate all independent cycles in the subgraph.
        try:
            cycles = nx.cycle_basis(subgraph)
        except Exception:
            cycles = []
        
        # Classify each cycle as a rigid ring or not.
        # Criteria for a rigid ring:
        #   - At least 5 atoms (i.e. not a 3- or 4-membered artefact).
        #   - Fewer than 3 rotatable bonds around the ring perimeter, meaning
        #     the ring is predominantly held together by non-rotatable bonds
        #     and should be treated as a single rigid unit.
        rigid_rings = []
        for cycle in cycles:
            if len(cycle) < 5:
                # Tiny cycles (3- or 4-membered) are not handled as rigid
                # sub-cores in the rotamer assignment.
                continue
            
            # Walk every consecutive pair of atoms in the cycle and count
            # how many of the ring bonds are rotatable.
            rotatable_bonds_in_cycle = 0
            for j in range(len(cycle)):
                node1 = cycle[j]
                node2 = cycle[(j + 1) % len(cycle)]  # wraps around at the last atom
                if self.has_edge(node1, node2) and self[node1][node2]['weight'] == 1:
                    rotatable_bonds_in_cycle += 1
            
            if rotatable_bonds_in_cycle < 3:
                rigid_rings.append(set(cycle))
        
        # Second pass: build a subgraph using only non-rotatable bonds and
        # run cycle detection again. This catches rigid rings whose perimeter
        # bonds are exclusively non-rotatable (e.g. aromatic rings), which
        # may have been missed or incorrectly excluded in the first pass due
        # to the rotatable-bond threshold.
        nonrot_subgraph = nx.Graph()
        nonrot_subgraph.add_nodes_from(branch_nodes)
        
        for node in branch_nodes:
            for neighbor in self.neighbors(node):
                if neighbor in branch_nodes and self[node][neighbor]['weight'] == 0:
                    nonrot_subgraph.add_edge(node, neighbor)
        
        try:
            nonrot_cycles = nx.cycle_basis(nonrot_subgraph)
            for cycle in nonrot_cycles:
                if len(cycle) >= 5:
                    cycle_set = set(cycle)
                    # Only add the cycle if it has not already been captured
                    # by the full-subgraph pass above.
                    if not any(cycle_set == ring for ring in rigid_rings):
                        rigid_rings.append(cycle_set)
        except Exception:
            pass
        
        return rigid_rings

    def _split_branch_at_rings(self, branch_nodes):
        """
        Split a branch into hierarchical groups only at true bifurcation points.

        A bifurcation occurs when a ring has 3+ rotatable bond connections
        (one from the core/base direction and 2+ extending outward).

        The algorithm is applied recursively:
          1. Find all bifurcation rings in *candidate_nodes*, sorted by
             distance to the current local core (closest first).
          2. At the closest bifurcation ring, measure the total number of
             rotatable bonds in each outward subbranch (all nodes reachable
             beyond the ring in that direction, going all the way to the ends).
          3. The dominant subbranch (most rotatable bonds) stays merged with
             the base segment and the ring → one group.
          4. Each minority subbranch becomes its own group (ring_nodes + all
             nodes of that subbranch, minus any bonds already in the dominant).
          5. The algorithm is then repeated on the dominant subbranch (with
             the bifurcation ring acting as the new local core) and on each
             minority subbranch independently.

        Parameters
        ----------
        branch_nodes : set[int]
            The nodes in the branch to split

        Returns
        -------
        hierarchical_groups : list[set[int]]
            List of node groups; each group carries a disjoint set of
            rotatable bonds.
        """
        # ------------------------------------------------------------------
        # Helpers
        # ------------------------------------------------------------------
        def collect_segment(start, candidate_nodes, exclude_nodes):
            """
            BFS over candidate_nodes from start, not crossing exclude_nodes.

            Parameters
            ----------
            start : int
                The starting node for the BFS
            candidate_nodes : set[int]
                The set of nodes to consider for the BFS
            exclude_nodes : set[int]
                The set of nodes to exclude from the BFS

            Returns
            -------
            segment : set[int]
                The set of nodes reached by the BFS
            """
            segment = set()
            stack = [start]
            visited = set()
            while stack:
                node = stack.pop()
                if node in visited or node in exclude_nodes or node not in candidate_nodes:
                    continue
                visited.add(node)
                segment.add(node)
                for nb in self.neighbors(node):
                    if nb not in visited and nb not in exclude_nodes and nb in candidate_nodes:
                        stack.append(nb)
            return segment

        def count_rot_bonds_in_segment(segment):
            """
            Count rotatable bonds whose both endpoints lie within segment.
            
            Parameters
            ----------
            segment : set[int]
                The set of nodes defining the segment

            Returns
            -------
            count : int
                The number of rotatable bonds within the segment
            """
            count = 0
            seen = set()
            for node in segment:
                for nb in self.neighbors(node):
                    if nb in segment and self[node][nb]['weight'] == 1:
                        key = frozenset([node, nb])
                        if key not in seen:
                            seen.add(key)
                            count += 1
            return count

        def dist_to_local_core(node, local_core_nodes):
            """
            Shortest unweighted path length from node to any local-core node.

            Parameters
            ----------
            node : int
                The node to measure distance from
            local_core_nodes : set[int]
                The nodes that act as the "core" for this recursion level
            Returns
            -------

            best : float
                The shortest unweighted path length from node to any local-core
                node, or infinity if no path exists
            """
            best = float('inf')
            for cn in local_core_nodes:
                try:
                    d = nx.shortest_path_length(self, node, cn)
                    if d < best:
                        best = d
                except nx.NetworkXNoPath:
                    pass
            return best

        # ------------------------------------------------------------------
        # Recursive splitter.
        #
        # candidate_nodes  – the set of branch nodes still to be assigned
        # local_core_nodes – nodes that act as the "core" for THIS recursion
        #                    level (originally the molecule's real core, then
        #                    the bifurcation ring found at the previous level)
        # base_accumulated – nodes accumulated on the core-side *above* this
        #                    recursion level (included in the dominant group)
        # ------------------------------------------------------------------
        def split_recursive(candidate_nodes, local_core_nodes, base_accumulated):
            """
            Returns a list of node-sets (groups).
            The first element is always the dominant chain from local_core
            toward the deepest dominant end.

            Parameters
            ----------
            candidate_nodes : set[int]
                The set of branch nodes still to be assigned
            local_core_nodes : set[int]
                Nodes that act as the "core" for THIS recursion level
                (originally the molecule's real core, then the bifurcation
                ring found at the previous level)
            base_accumulated : set[int]
                Nodes accumulated on the core-side *above* this recursion level
                (included in the dominant group)

            Returns
            -------
            groups : list[set[int]]
                A list of node-sets (groups). The first element is always the
                dominant chain from local_core toward the deepest dominant end.
            """
            # Detect rigid rings within the current candidate set and identify
            # which of them are bifurcation points, i.e. rings that have 3 or
            # more rotatable connections to the rest of the molecule (one
            # inward toward the core and at least two outward into distinct
            # sub-branches).
            rings = self._identify_rigid_rings(candidate_nodes)

            # --- Identify bifurcation rings ---------------------------------
            bifurcation_rings = []
            for ring_nodes in rings:
                # Collect every rotatable bond that exits the ring, either
                # toward the local-core or toward another candidate node that
                # lies outside the ring.
                rot_conns = []
                for rn in ring_nodes:
                    for nb in self.neighbors(rn):
                        is_local_core = nb in local_core_nodes
                        is_in_candidate_outside_ring = (
                            nb in candidate_nodes and nb not in ring_nodes)
                        if (is_local_core or is_in_candidate_outside_ring) \
                                and self[rn][nb]['weight'] == 1:
                            rot_conns.append((rn, nb))
                if len(rot_conns) >= 3:
                    bifurcation_rings.append((ring_nodes, rot_conns))

            if not bifurcation_rings:
                # No bifurcation ring found: the entire candidate set forms a
                # single linear (or branching-but-not-ring-split) group.
                group = base_accumulated | candidate_nodes
                return [group]

            # --- Sort bifurcation rings by distance to local core -----------
            def ring_distance(ring_and_conns):
                rn_set, _ = ring_and_conns
                return min(dist_to_local_core(n, local_core_nodes)
                           for n in rn_set)

            bifurcation_rings.sort(key=ring_distance)
            nearest_ring_nodes, connections = bifurcation_rings[0]

            # --- Identify the core-direction connection ----------------------
            # The inward connection is the rotatable bond that leads back
            # toward the local core. Prefer a direct connection (the external
            # node is already in local_core_nodes); fall back to the
            # connection whose external node is closest to the local core.
            inward_conn = None
            for rn, nb in connections:
                if nb in local_core_nodes:
                    inward_conn = (rn, nb)
                    break
            if inward_conn is None:
                min_d = float('inf')
                for rn, nb in connections:
                    d = dist_to_local_core(nb, local_core_nodes)
                    if d < min_d:
                        min_d = d
                        inward_conn = (rn, nb)

            # --- Base segment: nodes on the core-side of the ring -----------
            # Collect all candidate nodes reachable from the inward external
            # node without crossing the ring. These nodes form the "stem"
            # that connects the local core to the bifurcation ring and will
            # be merged into the dominant group.
            inward_external = inward_conn[1]
            if inward_external in candidate_nodes:
                base_segment = collect_segment(
                    inward_external,
                    candidate_nodes=candidate_nodes,
                    exclude_nodes=nearest_ring_nodes)
            else:
                # The ring is directly attached to the local core with no
                # intermediate nodes.
                base_segment = set()

            # --- Outward connections (all except the inward one) ------------
            # These are the rotatable bonds that exit the ring toward distinct
            # sub-branches (i.e. the "forks" that justify the split).
            outward_conns = [
                (rn, nb) for rn, nb in connections
                if nb not in local_core_nodes
                and (rn, nb) != inward_conn
            ]

            # --- Collect full subbranch for each outward connection ---------
            # A "full subbranch" is the set of all candidate nodes reachable
            # from an outward exit node without re-crossing the bifurcation
            # ring or the base segment. The total number of rotatable bonds
            # within each subbranch (including the exit bond itself) is used
            # to rank them.
            excluded_base = nearest_ring_nodes | base_segment
            outward_subbranches = []
            for rn, nb in outward_conns:
                if nb not in candidate_nodes:
                    continue
                full_sub = collect_segment(
                    nb,
                    candidate_nodes=candidate_nodes,
                    exclude_nodes=excluded_base)
                rot_count = count_rot_bonds_in_segment(full_sub)
                # Include the rotatable bond from the ring into this subbranch
                # in the count so that single-bond branches are not penalised.
                if self[rn][nb]['weight'] == 1:
                    rot_count += 1
                outward_subbranches.append((nb, full_sub, rot_count))

            if not outward_subbranches:
                # All outward connections led outside the candidate set, so
                # there is nothing to split – treat everything as one group.
                group = base_accumulated | base_segment | nearest_ring_nodes | candidate_nodes
                return [group]

            # --- Select dominant subbranch (most rotatable bonds) -----------
            # The dominant subbranch is the one with the highest total number
            # of rotatable bonds; it stays merged with the core-side segment
            # and the ring. All other subbranches become separate groups.
            dom_idx = max(range(len(outward_subbranches)),
                          key=lambda i: outward_subbranches[i][2])
            dom_start, dom_sub, dom_rot = outward_subbranches[dom_idx]

            # --- Build dominant group & recurse on it -----------------------
            # The new "local core" for the dominant recursion is the ring itself.
            dominant_base = base_accumulated | base_segment | nearest_ring_nodes
            dominant_groups = split_recursive(
                candidate_nodes=dom_sub,
                local_core_nodes=nearest_ring_nodes,
                base_accumulated=dominant_base)

            # The first dominant_group absorbs dominant_base (already done
            # inside the recursive call via base_accumulated).

            # --- Build minority groups & recurse on each --------------------
            # Each minority subbranch becomes its own group (or further splits
            # if it itself contains bifurcation rings). The bifurcation ring is
            # passed as the local-core for these sub-recursions so that bond
            # ordering inside each minority branch is measured from the ring
            # outward.
            minority_groups = []
            for i, (sub_start, sub_nodes, sub_rot) in enumerate(outward_subbranches):
                if i == dom_idx:
                    continue
                min_groups = split_recursive(
                    candidate_nodes=sub_nodes,
                    local_core_nodes=nearest_ring_nodes,
                    base_accumulated=nearest_ring_nodes.copy())
                minority_groups.extend(min_groups)

            all_groups = dominant_groups + minority_groups
            return all_groups

        # ------------------------------------------------------------------
        # Entry point
        # ------------------------------------------------------------------
        # Find the node in branch_nodes that is directly connected to core
        local_core = set(self.core_nodes)

        # Start the recursion
        groups = split_recursive(
            candidate_nodes=branch_nodes,
            local_core_nodes=local_core,
            base_accumulated=set())

        return groups

    def get_rotamers(self):
        """
        It builds the RotamerLibrary object with hierarchical branch splitting.

        Returns
        -------
        rotamers : list[list]
            The list of rotamers grouped by the branch they belong to
        """
        resolution = self.molecule.rotamer_resolution

        assert len(self.core_nodes) > 0, 'No core nodes were found'

        # Remove core nodes to isolate the peripheral branches. Each connected
        # component that remains is an independent branch originating from the
        # core.
        branch_graph = deepcopy(self)

        for node in self.core_nodes:
            branch_graph.remove_node(node)

        branch_groups = list(nx.connected_components(branch_graph))
        
        # Attempt to split each branch at bifurcating rings. When a branch
        # contains no ring-based bifurcations, _split_branch_at_rings returns
        # the original node set unchanged (wrapped in a single-element list),
        # so the extend call is always safe.
        hierarchical_branch_groups = []
        for branch_nodes in branch_groups:
            split_groups = self._split_branch_at_rings(branch_nodes)
            hierarchical_branch_groups.extend(split_groups)

        # Map every rotatable bond to its branch group, then determine the
        # core-adjacent atom for each group so bonds can be sorted from the
        # core outward.
        rot_bonds_per_group = self._get_rot_bonds_per_group(hierarchical_branch_groups)

        core_atom_per_group = self._get_core_atom_per_group(rot_bonds_per_group)

        # Sort rotatable bonds within each group by increasing graph distance
        # from the group's core atom, ensuring a core-to-tip ordering.
        distances = dict(nx.shortest_path_length(self))
        sorted_rot_bonds_per_group = self._get_sorted_bonds_per_group(
            core_atom_per_group, rot_bonds_per_group, distances)

        rotamers = list()

        # Build Rotamer objects from the ordered bond lists; skip any group
        # that ended up with no rotatable bonds after filtering.
        for rot_bonds in sorted_rot_bonds_per_group:
            branch_rotamers = list()
            for (atom1_index, atom2_index) in rot_bonds:
                rotamer = Rotamer(atom1_index, atom2_index, resolution)
                branch_rotamers.append(rotamer)

            if len(branch_rotamers) > 0:
                rotamers.append(branch_rotamers)

        # Place the group with the most rotatable bonds first so that PELE
        # uses the longest chain as the primary branch.
        rotamers.sort(key=len, reverse=True)

        return rotamers

    @property
    def molecule(self):
        """
        The peleffy's Molecule.

        Returns
        -------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object
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


class MolecularGraphWithConstrainedCore(MolecularGraph):
    """
    It represents the structure of a Molecule as a networkx.Graph.
    """

    def __init__(self, molecule, atom_constraints):
        """
        It initializes a MolecularGraph object.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        atom_constraint : list[int or str]
            It defines the list of atoms to constrain in the core, thus,
            the core will be forced to contain them. Atoms can be specified
            through integers that match the atom index or strings that
            match with the atom PDB name

        Raises
        ------
        ValueError
            If the supplied array of atom constraints is empty
            If the PDB atom name in atom_constraint does not match with
            any atom in the molecule
        TypeError
            If the atom_constraint is of invalid type
        """
        if len(atom_constraints) == 0:
            raise ValueError('Supplied empty array of atom constraints')

        self._constraint_indices = list()

        for atom_constraint in atom_constraints:
            if isinstance(atom_constraint, int):
                self._constraint_indices.append(atom_constraint)
            elif isinstance(atom_constraint, str):
                atom_names = molecule.get_pdb_atom_names()
                for index, name in enumerate(atom_names):
                    if name == atom_constraint:
                        self._constraint_indices.append(index)
                        break
                else:
                    raise ValueError('Supplied PDB atom name '
                                     + '\'{}\''.format(atom_constraint)
                                     + 'is missing in molecule')
            else:
                raise TypeError('Invalid type for the atom_constraint')

        super().__init__(molecule)

        self._safety_check()

    def _build_core_nodes(self):
        """
        It builds the list of core nodes
        """

        from networkx.algorithms.shortest_paths.generic import \
            shortest_path_length

        # Force core to contain constrained atoms
        core_nodes = [index for index in self.constraint_indices]

        # Calculate graph distances according to weight values
        weighted_distances = dict(shortest_path_length(self, weight="weight"))

        # Add also all atoms at 0 distance with respect to constrained
        # atom into the core
        for node in self.nodes:
            for constraint_index in self.constraint_indices:
                d = weighted_distances[constraint_index][node]
                if d == 0 and node not in core_nodes:
                    core_nodes.append(node)

        self._core_nodes = core_nodes

    def _safety_check(self):
        """Perform a safety check on the atom constraints."""
        if len(self.constraint_indices) < 2:
            return

        safe_indices = set()
        for i, cidx1 in enumerate(self.constraint_indices):
            if cidx1 in safe_indices:
                continue
            for cidx2 in self.constraint_indices:
                if cidx2 in self.neighbors(cidx1):
                    safe_indices.add(cidx1)
                    safe_indices.add(cidx2)
                    break
            else:
                raise ValueError('All atoms in atom constraints must be '
                                 + 'adjacent and atom '
                                 + self.constraint_names[i].strip()
                                 + ' is not')

    @property
    def constraint_indices(self):
        """
        The indices of atoms to constraint to the core.

        Returns
        -------
        constraint_indices : list[int]
            List of atom indices
        """
        return self._constraint_indices

    @property
    def constraint_names(self):
        """
        The names of atoms to constraint to the core.

        Returns
        -------
        constraint_names : list[str]
            List of atom names
        """
        atom_names = self.molecule.get_pdb_atom_names()

        constraint_names = list()
        for index in self.constraint_indices:
            constraint_names.append(atom_names[index])

        return constraint_names


class CoreLessMolecularGraph(MolecularGraph):
    """
    It represents the structure of a Molecule as a networkx.Graph with
    an empty atom core. This molecular graph is only valid for
    alchemy applications.
    """

    def __init__(self, molecule):
        """
        It initializes a CoreLessMolecularGraph object.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        nx.Graph.__init__(self)
        self._molecule = molecule
        self._compute_rotamer_graph()
        self._core_nodes = []

    def get_rotamers(self):
        """
        It builds the RotamerLibrary object with hierarchical branch splitting.

        Returns
        -------
        rotamers : list[list]
            The list of rotamers grouped by the branch they belong to
        """
        resolution = self.molecule.rotamer_resolution

        branch_graph = deepcopy(self)

        branch_groups = list(nx.connected_components(branch_graph))
        
        # Apply hierarchical splitting to each branch
        hierarchical_branch_groups = []
        for branch_nodes in branch_groups:
            split_groups = self._split_branch_at_rings(branch_nodes)
            hierarchical_branch_groups.extend(split_groups)

        rot_bonds_per_group = self._get_rot_bonds_per_group(hierarchical_branch_groups)

        rotamers = list()

        for group_id, rot_bonds in enumerate(rot_bonds_per_group):
            branch_rotamers = list()
            for (atom1_index, atom2_index) in rot_bonds:
                rotamer = Rotamer(atom1_index, atom2_index, resolution)
                branch_rotamers.append(rotamer)

            if len(branch_rotamers) > 0:
                rotamers.append(branch_rotamers)

        return rotamers
