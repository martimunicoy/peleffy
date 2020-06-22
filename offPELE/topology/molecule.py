
# Global imports
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

from .rotamer import RotamerLibrary, Rotamer
from offPELE.utils.toolkits import (AmberToolkitWrapper,
                                    RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper)


class Atom(object):
    def __init__(self, index=-1, core=None, OPLS_type=None, PDB_name=None,
                 unknown=None, z_matrix_x=None, z_matrix_y=None,
                 z_matrix_z=None, sigma=None, epsilon=None, charge=None,
                 born_radius=None, SASA_radius=None, nonpolar_gamma=None,
                 nonpolar_alpha=None):
        self.index = index
        self.core = core
        self.OPLS_type = OPLS_type
        self.PDB_name = PDB_name
        self.unknown = unknown
        self.z_matrix_x = z_matrix_x
        self.z_matrix_y = z_matrix_y
        self.z_matrix_z = z_matrix_z
        self.sigma = sigma
        self.epsilon = epsilon
        self.charge = charge
        self.born_radius = born_radius  # Rad. Non Polar SGB
        self.SASA_radius = SASA_radius  # Rad. Non Polar Type
        self.nonpolar_gamma = nonpolar_gamma  # SGB Non Polar gamma
        self.nonpolar_alpha = nonpolar_alpha  # SGB Non Polar type

    def set_as_core(self):
        self.core = True

    def set_as_branch(self):
        self.core = False


class Bond(object):
    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 spring_constant=None, eq_dist=None):
        self.index = index
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.spring_constant = spring_constant
        self.eq_dist = eq_dist


class Molecule(object):
    def __init__(self, path=None):
        if isinstance(path, str):
            from pathlib import Path
            extension = Path(path).suffix
            extension = extension.strip('.')
            if extension == 'pdb':
                self._initialize_from_pdb(path)
            else:
                raise ValueError(
                    '{} is not a valid extension'.format(extension))
        else:
            self._initialize()

    def _initialize(self):
        self._name = ''
        self._forcefield = ''
        self._atoms = list()
        self._bonds = list()
        self._angles = list()
        self._dihedrals = list()
        self._impropers = list()
        self._rdkit_molecule = None
        self._off_molecule = None
        self._rotamer_library = None

    def _initialize_from_pdb(self, path):
        self._initialize()
        print(' - Loading molecule from RDKit')

        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_pdb(path)

        # RDKit must generate stereochemistry specifically from 3D coords
        rdkit_toolkit.assign_stereochemistry_from_3D(self)

        openforcefield_toolkit = OpenForceFieldToolkitWrapper()

        self._off_molecule = openforcefield_toolkit.from_rdkit(self)

        name = Path(path).stem
        if len(name) > 2:
            self.set_name(name)

    def set_name(self, name):
        if isinstance(name, str) and len(name) > 2:
            name = name[0:3].upper()
            self._name = name

            if self.off_molecule:
                self.off_molecule.name = name

    def parameterize(self, forcefield):
        if not self.off_molecule and not self.rdkit_molecule:
            raise Exception('OpenForceField molecule was not initialized '
                            + 'correctly')

        print(' - Loading forcefield')
        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        parameters = openforcefield_toolkit.get_parameters_from_forcefield(
            forcefield, self)

        self.parameters = parameters

        print(' - Computing partial charges with am1bcc')
        self._calculate_am1bcc_charges()

        self._build_atoms()

        self._build_bonds()

        self._build_angles()

        self._build_dihedrals()

        self._build_impropers()

    def _compute_rotamer_graph(self, rot_bonds_atom_ids):
        import networkx as nx
        graph = nx.Graph()

        for atom in self.rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            graph.add_node(atom.GetIdx(), pdb_name=pdb_info.GetName(),
                           nrot_neighbors=list())

        for bond in self.rdkit_molecule.GetBonds():
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()
            if ((atom1, atom2) in rot_bonds_atom_ids
                    or (atom2, atom2) in rot_bonds_atom_ids):
                rotatable = True
            else:
                rotatable = False
                graph.nodes[atom1]['nrot_neighbors'].append(atom2)
                graph.nodes[atom2]['nrot_neighbors'].append(atom1)

            graph.add_edge(bond.GetBeginAtomIdx(),
                           bond.GetEndAtomIdx(),
                           weight=int(rotatable))

        for i, j in rot_bonds_atom_ids:
            graph[i][j]['weight'] = 1
            graph.nodes[i]['rotatable'] = True
            graph.nodes[j]['rotatable'] = True

        return graph

    def _set_core(self, graph, rot_bonds_atom_ids):
        def get_all_nrot_neighbors(atom_id, visited_neighbors):
            if atom_id in visited_neighbors:
                return visited_neighbors
            visited_neighbors.add(atom_id)
            nrot_neighbors = graph.nodes[atom_id]['nrot_neighbors']
            for nrot_neighbor in nrot_neighbors:
                visited_neighbors = get_all_nrot_neighbors(
                    nrot_neighbor, visited_neighbors)
            return visited_neighbors

        from networkx.algorithms.shortest_paths.generic import \
            shortest_path_length
        from networkx.algorithms.distance_measures import eccentricity

        # Calculate graph distances according to weight values
        weighted_distances = dict(shortest_path_length(graph, weight="weight"))

        # Calculate eccentricites using weighted distances
        eccentricities = eccentricity(graph, sp=weighted_distances)

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
            centered_node_groups.append(get_all_nrot_neighbors(node, set()))

        # In case of more than one group, core will be the largest
        node_group = sorted(centered_node_groups, key=len, reverse=True)[0]

        for atom in self.atoms:
            if atom.index in node_group:
                atom.set_as_core()
            else:
                atom.set_as_branch()

    def _get_rot_bonds_per_group(self, graph, branch_groups):
        rot_bonds_per_group = list()
        for group in branch_groups:
            rot_bonds = list()
            visited_bonds = set()
            for node in group:
                bonds = graph.edges(node)
                for bond in bonds:
                    if bond in visited_bonds:
                        continue
                    if graph[bond[0]][bond[1]]['weight'] == 1:
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

    def _build_rotamer_library(self, graph, resolution):
        import networkx as nx

        core_atoms = set()
        for atom in self.atoms:
            if atom.core:
                core_atoms.add(atom)
        core_indexes = [atom.index for atom in core_atoms]

        assert len(core_atoms) > 0, 'No core atoms were found'

        branch_graph = deepcopy(graph)

        for core_atom in core_atoms:
            branch_graph.remove_node(core_atom.index)

        branch_groups = list(nx.connected_components(branch_graph))

        rot_bonds_per_group = self._get_rot_bonds_per_group(
            graph, branch_groups)

        core_atom_per_group = self._get_core_atom_per_group(
            rot_bonds_per_group, core_indexes)

        distances = dict(nx.shortest_path_length(graph))

        sorted_rot_bonds_per_group = self._get_sorted_bonds_per_group(
            core_atom_per_group, rot_bonds_per_group, distances)

        self._rotamer_library = RotamerLibrary(self.name)

        # PELE needs underscores instead of whitespaces
        pdb_atom_names = [name.replace(' ', '_',)
                          for name in self.get_pdb_atom_names()]

        for group_id, rot_bonds in enumerate(sorted_rot_bonds_per_group):
            for (atom1_index, atom2_index) in rot_bonds:
                atom1_name = pdb_atom_names[atom1_index]
                atom2_name = pdb_atom_names[atom2_index]
                rotamer = Rotamer(atom1_name, atom2_name, resolution)
                self.rotamer_library.add_rotamer(rotamer, group_id)

    def build_rotamer_library(self, resolution):
        self._assert_parameterized()

        try:
            from rdkit import Chem
        except ImportError:
            raise Exception('RDKit Python API not found')

        print(' - Generating rotamer library')

        # Fins rotatable bond ids as in Lipinski module in RDKit
        # https://github.com/rdkit/rdkit/blob/1bf6ef3d65f5c7b06b56862b3fb9116a3839b229/rdkit/Chem/Lipinski.py#L47
        rot_bonds_atom_ids = self._rdkit_molecule.GetSubstructMatches(
            Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]'))

        graph = self._compute_rotamer_graph(rot_bonds_atom_ids)

        self._set_core(graph, rot_bonds_atom_ids)

        self._build_rotamer_library(graph, resolution)

    def plot_rotamer_graph(self):
        self._assert_parameterized()

        try:
            from rdkit import Chem
        except ImportError:
            raise Exception('RDKit Python API not found')

        # Fins rotatable bond ids as in Lipinski module in RDKit
        # https://github.com/rdkit/rdkit/blob/1bf6ef3d65f5c7b06b56862b3fb9116a3839b229/rdkit/Chem/Lipinski.py#L47
        rot_bonds_atom_ids = self._rdkit_molecule.GetSubstructMatches(
            Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]'))

        graph = self._compute_rotamer_graph(rot_bonds_atom_ids)

        rot_edges = [(u, v) for (u, v, d) in graph.edges(data=True)
                     if d['weight'] == 1]
        nrot_edges = [(u, v) for (u, v, d) in graph.edges(data=True)
                      if d['weight'] == 0]

        import networkx as nx

        pos = nx.circular_layout(graph)

        nx.draw_networkx_nodes(graph, pos, node_size=400)
        nx.draw_networkx_edges(graph, pos, edgelist=rot_edges,
                               width=4)
        nx.draw_networkx_edges(graph, pos, edgelist=nrot_edges,
                               width=4, alpha=0.5, edge_color='b',
                               style='dashed')
        nx.draw_networkx_labels(graph, pos, font_size=10,
                                font_family='sans-serif')

        import matplotlib.pyplot as plt
        plt.axis('off')
        plt.show()

    def _assert_parameterized(self):
        try:
            assert self.off_molecule is not None
        except AssertionError:
            raise Exception('Molecule not parameterized')

    def _calculate_am1bcc_charges(self):
        ambertoolkit = AmberToolkitWrapper()

        charges = ambertoolkit.compute_partial_charges_am1bcc(self)

        self.off_molecule.partial_charges = charges

    def _build_atoms(self):
        # PELE needs underscores instead of whitespaces
        pdb_atom_names = [name.replace(' ', '_',)
                          for name in self.get_pdb_atom_names()]

        for index, atom in enumerate(self.off_molecule.atoms):
            atom = Atom(index=index, PDB_name=pdb_atom_names[index],
                        OPLS_type=None, unknown=None, z_matrix_x=None,
                        z_matrix_y=None, z_matrix_z=None,
                        sigma=None,
                        epsilon=None,
                        charge=self.off_molecule.partial_charges[index],
                        born_radius=None,
                        SASA_radius=None,
                        nonpolar_gamma=None,
                        nonpolar_alpha=None)
            self._add_atom(atom)

    def _add_atom(self, atom):
        self._atoms.append(atom)

    def _build_bonds(self):
        pass

    def _build_angles(self):
        pass

    def _build_dihedrals(self):
        pass

    def _build_impropers(self):
        pass

    def get_pdb_atom_names(self):
        self._assert_parameterized()

        pdb_atom_names = list()

        for atom in self.rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            pdb_atom_names.append(pdb_info.GetName())

        return pdb_atom_names

    def to_impact(self, path):
        pass

    @property
    def off_molecule(self):
        return self._off_molecule

    @property
    def rdkit_molecule(self):
        return self._rdkit_molecule

    @property
    def rotamer_library(self):
        return self._rotamer_library

    @property
    def name(self):
        return self._name

    @property
    def atoms(self):
        return self._atoms
