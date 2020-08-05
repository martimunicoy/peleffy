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

    def __init__(self, atom1, atom2, resolution=30):
        """
        It initiates an Rotamer object.

        Parameters
        ----------
        atom1 : str
            The name of the first atom involved in the rotamer
        atom2 : str
            The name of the second atom involved in the rotamer
        resolution : float
            The resolution to discretize the rotamer's conformational space
        """
        self._atom1 = atom1
        self._atom2 = atom2
        self._resolution = resolution

    @property
    def atom1(self):
        """
        Rotamer's atom1 name.

        Returns
        -------
        atom1 : str
            The name of the first atom involved in this Rotamer object
        """
        return self._atom1

    @property
    def atom2(self):
        """
        Rotamer's atom2 name.

        Returns
        -------
        atom2 : str
            The name of the second atom involved in this Rotamer object
        """
        return self._atom2

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

    def __init__(self, residue_name='LIG'):
        """
        It initiates a RotamerLibrary object.

        Parameters
        ----------
        residue_name : str
            The name of the residue this RotamerLibrary belongs to
        """
        self._residue_name = residue_name
        self._rotamers = defaultdict(list)

    def add_rotamer(self, rotamer, group_id):
        """
        It adds a rotamer to this RotamerLibrary.

        Parameters
        ----------
        rotamer : a Rotamer
            The Rotamer object to add to this RotamerLibrary
        group_id : int
            The index of the group where the rotamer belongs to
        """
        self._rotamers[group_id].append(rotamer)

    def to_file(self, path):
        """
        It writes this RotamerLibrary to a file.

        Parameters
        ----------
        path : str
            Path to save the RotamerLibrary to
        """
        with open(path, 'w') as file:
            file.write('rot assign res {} &\n'.format(self.residue_name))
            for i, group in enumerate(self.rotamers.keys()):
                if i > 0:
                    file.write('     newgrp &\n')
                for rotamer in self.rotamers[group]:
                    file.write('   sidelib FREE{} {} {} &\n'.format(
                        rotamer.resolution, rotamer.atom1, rotamer.atom2))

    @property
    def residue_name(self):
        """
        The name of the RotamerLibrary's residue.

        Returns
        -------
        residue_name : str
            The name of the residue of this RotamberLibrary object.
        """
        return self._residue_name

    @property
    def rotamers(self):
        """
        RotamerLibrary's rotamers.

        Returns
        -------
        rotamers : dict[int, Rotamer]
            The rotamers of this RotamerLibrary object, grouped by the
            group they belong to
        """
        return self._rotamers


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
        self._compute_rotamer_graph(molecule)

    def _compute_rotamer_graph(self, molecule):
        """
        It initializes the netwrokx.Graph with a Molecule object.

        Parameters
        ----------
        molecule : An offpele.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        rdkit_toolkit = RDKitToolkitWrapper()
        rot_bonds_atom_ids = \
            rdkit_toolkit.get_atom_ids_with_rotatable_bonds(molecule)

        rdkit_molecule = molecule.rdkit_molecule

        for atom in rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            self.add_node(atom.GetIdx(), pdb_name=pdb_info.GetName(),
                          nrot_neighbors=list())

        for bond in rdkit_molecule.GetBonds():
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()
            if ((atom1, atom2) in rot_bonds_atom_ids
                    or (atom2, atom2) in rot_bonds_atom_ids):
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

    def set_core(self):
        """
        It sets the core of the molecule to minimize the amount of consecutive
        rotamers as much as possible.
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
        already_visited = set()
        centered_node_groups = list()
        for node in centered_nodes:
            if node in already_visited:
                continue
            centered_node_groups.append(get_all_nrot_neighbors(self, node,
                                                               set()))

        # In case of more than one group, core will be the largest
        node_group = sorted(centered_node_groups, key=len, reverse=True)[0]

        for atom in self.molecule.atoms:
            if atom.index in node_group:
                atom.set_as_core()
            else:
                atom.set_as_branch()

    def set_parents(self):
        """
        It sets the parent of each atom according to the molecular graph.
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

    def _get_core_atom_per_group(self, rot_bonds_per_group, core_indexes):
        """
        It obtains the core atom for each group.

        Parameters
        ----------
        rot_bonds_per_group : list[tuple[int, int]]
            The atom ids of all the graph's edges that belong to a rotatable
            bond
        core_indexes : list[int]
            The atom ids of atoms in the core

        Returns
        -------
        core_atom_per_group : list[int]
            The atom id of the atom that belongs to the core for each branch
        """
        core_atom_per_group = list()
        for rot_bonds in rot_bonds_per_group:
            for (a1, a2) in rot_bonds:
                if a1 in core_indexes:
                    core_atom_per_group.append(a1)
                    break
                elif a2 in core_indexes:
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
                                         n_rot_bonds_to_ignore):
        """
        It ignores a certain number of terminal rotatable bonds of each
        group.

        Parameters
        ----------
        sorted_rot_bonds_per_group : list[list]
            The rotatable bonds per group, sorted in increasing order by
            their distance with respect to the corresponding core atom
        n_rot_bonds_to_ignore : int
            The number of terminal rotatable bonds to ignore in each group

        Returns
        -------
        filtered_rot_bonds_per_group : list[list]
            The filtered rotatable bonds per group that are obtained
        """
        if n_rot_bonds_to_ignore == 0:
            return sorted_rot_bonds_per_group

        filtered_rot_bonds_per_group = list()

        for rot_bonds in sorted_rot_bonds_per_group:
            filtered_rot_bonds_per_group.append(
                rot_bonds[:-n_rot_bonds_to_ignore])

        return filtered_rot_bonds_per_group

    def build_rotamer_library(self, resolution=30, n_rot_bonds_to_ignore=1):
        """
        It builds the RotamerLibrary object.

        Parameters
        ----------
        resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        n_rot_bonds_to_ignore : int
            The number of terminal rotatable bonds to ignore when
            building the rotamer library. Default is 1

        Returns
        -------
        rotamer_library : a RotamerLibrary object
            The RotamerLibrary for the supplied Molecule object.
        """
        core_atoms = set()
        for atom in self.molecule.atoms:
            if atom.core:
                core_atoms.add(atom)
        core_indexes = [atom.index for atom in core_atoms]

        assert len(core_atoms) > 0, 'No core atoms were found'

        branch_graph = deepcopy(self)

        for core_atom in core_atoms:
            branch_graph.remove_node(core_atom.index)

        branch_groups = list(nx.connected_components(branch_graph))

        rot_bonds_per_group = self._get_rot_bonds_per_group(branch_groups)

        core_atom_per_group = self._get_core_atom_per_group(
            rot_bonds_per_group, core_indexes)

        distances = dict(nx.shortest_path_length(self))

        sorted_rot_bonds_per_group = self._get_sorted_bonds_per_group(
            core_atom_per_group, rot_bonds_per_group, distances)

        filtered_rot_bonds_per_group = self._ignore_terminal_rotatable_bonds(
            sorted_rot_bonds_per_group, n_rot_bonds_to_ignore)

        rotamer_library = RotamerLibrary(self.molecule.name)

        # PELE needs underscores instead of whitespaces
        pdb_atom_names = [name.replace(' ', '_',)
                          for name in self.molecule.get_pdb_atom_names()]

        for group_id, rot_bonds in enumerate(filtered_rot_bonds_per_group):
            for (atom1_index, atom2_index) in rot_bonds:
                atom1_name = pdb_atom_names[atom1_index]
                atom2_name = pdb_atom_names[atom2_index]
                rotamer = Rotamer(atom1_name, atom2_name, resolution)
                rotamer_library.add_rotamer(rotamer, group_id)

        return rotamer_library

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
