# Global imports
from collections import defaultdict
import networkx as nx
from copy import deepcopy

from offPELE.utils.toolkits import RDKitToolkitWrapper


class Rotamer(object):
    def __init__(self, atom1, atom2, resolution=30):
        self._atom1 = atom1
        self._atom2 = atom2
        self._resolution = resolution

    @property
    def atom1(self):
        return self._atom1

    @property
    def atom2(self):
        return self._atom2

    @property
    def resolution(self):
        return self._resolution


class RotamerLibrary(object):
    def __init__(self, residue_name='LIG'):
        self._residue_name = residue_name
        self._rotamers = defaultdict(list)

    def add_rotamer(self, rotamer, group_id):
        self._rotamers[group_id].append(rotamer)

    def to_file(self, path):
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
        return self._residue_name

    @property
    def rotamers(self):
        return self._rotamers


class MolecularGraph(nx.Graph):
    def __init__(self, molecule):
        super().__init__(self)
        self._molecule = molecule
        self._compute_rotamer_graph(molecule)

    def _compute_rotamer_graph(self, molecule):
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
        def get_all_nrot_neighbors(self, atom_id, visited_neighbors):
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
        def recursive_child_visitor(parent, already_visited=set()):
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

        # Assert absolut parent is the only with None parent value
        assert parent.parent is None and \
            sum([int(a.parent is not None) for a in self.molecule.atoms]) \
            == len(self.molecule.atoms) - 1, 'Found descendant without parent'

    def _get_rot_bonds_per_group(self, branch_groups):
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

    def build_rotamer_library(self, resolution):
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

        rotamer_library = RotamerLibrary(self.molecule.name)

        # PELE needs underscores instead of whitespaces
        pdb_atom_names = [name.replace(' ', '_',)
                          for name in self.molecule.get_pdb_atom_names()]

        for group_id, rot_bonds in enumerate(sorted_rot_bonds_per_group):
            for (atom1_index, atom2_index) in rot_bonds:
                atom1_name = pdb_atom_names[atom1_index]
                atom2_name = pdb_atom_names[atom2_index]
                rotamer = Rotamer(atom1_name, atom2_name, resolution)
                rotamer_library.add_rotamer(rotamer, group_id)

        return rotamer_library

    @property
    def molecule(self):
        return self._molecule
