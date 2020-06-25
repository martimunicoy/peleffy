
# Global imports
from collections import defaultdict
from copy import deepcopy
from pathlib import Path

from .rotamer import RotamerLibrary, Rotamer
from .topology import Bond, Angle, OFFProper, OFFImproper
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
        self._forcefield = None
        self._atoms = list()
        self._bonds = list()
        self._angles = list()
        self._propers = list()
        self._OFF_propers = list()
        self._impropers = list()
        self._OFF_impropers = list()
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
        # TODO Is there a way to retrieve the name of the OFF's ForceField object?
        if isinstance(forcefield, str):
            self._forcefield = Path(forcefield).stem

        print(' - Computing partial charges with am1bcc')
        self._calculate_am1bcc_charges()

        self._build_atoms()

        self._build_bonds()

        self._build_angles()

        self._build_propers()

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
        amber_toolkit = AmberToolkitWrapper()

        charges = amber_toolkit.compute_partial_charges_am1bcc(self)

        self.off_molecule.partial_charges = charges

    def _build_atoms(self):
        # PELE needs underscores instead of whitespaces
        pdb_atom_names = {(i, ): name.replace(' ', '_',)
                          for i, name in enumerate(self.get_pdb_atom_names())}

        # TODO should we assign an OPLS type? How can we do this with OFF?
        OPLS_types = {i: None
                      for i in self.parameters.get_vdW_parameters().keys()}

        # TODO Which is the purpose of unknown value? Is it important?
        unknowns = {i: None
                    for i in self.parameters.get_vdW_parameters().keys()}

        # TODO Create z-matrix from 3D coordinates
        z_matrix_xs, z_matrix_ys, z_matrix_zs = (
            {i: None for i in self.parameters.get_vdW_parameters().keys()},
            {i: None for i in self.parameters.get_vdW_parameters().keys()},
            {i: None for i in self.parameters.get_vdW_parameters().keys()})

        sigmas = self.parameters.get_vdW_sigmas()

        if all([sigma is None for sigma in sigmas.values()]):
            sigmas = self.parameters.get_vdW_sigmas_from_rmin_halves()

        epsilons = self.parameters.get_vdW_epsilons()

        # TODO Find a way to assign implicit solvent parameters to atoms with OFF
        born_radii = {i: None
                      for i in self.parameters.get_vdW_parameters().keys()}
        SASA_radii = {i: None
                      for i in self.parameters.get_vdW_parameters().keys()}
        nonpolar_gammas = {i: None for i in
                           self.parameters.get_vdW_parameters().keys()}
        nonpolar_alphas = {i: None for i in
                           self.parameters.get_vdW_parameters().keys()}

        for index in self.parameters.get_vdW_parameters().keys():
            atom = Atom(index=index, PDB_name=pdb_atom_names[index],
                        OPLS_type=OPLS_types[index],
                        unknown=unknowns[index],
                        z_matrix_x=z_matrix_xs[index],
                        z_matrix_y=z_matrix_ys[index],
                        z_matrix_z=z_matrix_zs[index],
                        sigma=sigmas[index],
                        epsilon=epsilons[index],
                        charge=self.off_molecule.partial_charges[index],
                        born_radius=born_radii[index],
                        SASA_radius=SASA_radii[index],
                        nonpolar_gamma=nonpolar_gammas[index],
                        nonpolar_alpha=nonpolar_alphas[index])
            self._add_atom(atom)

    def _add_atom(self, atom):
        self._atoms.append(atom)

    def _build_bonds(self):
        lengths = self.parameters.get_bond_lengths()

        ks = self.parameters.get_bond_ks()

        bond_indexes = self.parameters.get_bond_parameters().keys()

        for index, atom_indexes in enumerate(bond_indexes):
            (atom1_idx, atom2_idx) = atom_indexes
            bond = Bond(index=index, atom1_idx=atom1_idx, atom2_idx=atom2_idx,
                        spring_constant=ks[atom_indexes],
                        eq_dist=lengths[atom_indexes])
            self._add_bond(bond)

    def _add_bond(self, bond):
        self._bonds.append(bond)

    def _build_angles(self):
        angles = self.parameters.get_angle_angles()

        ks = self.parameters.get_angle_ks()

        angle_indexes = self.parameters.get_angle_parameters().keys()

        for index, atom_indexes in enumerate(angle_indexes):
            atom1_idx, atom2_idx, atom3_idx = atom_indexes
            angle = Angle(index=index, atom1_idx=atom1_idx,
                          atom2_idx=atom2_idx, atom3_idx=atom3_idx,
                          spring_constant=ks[atom_indexes],
                          eq_angle=angles[atom_indexes])
            self._add_angle(angle)

    def _add_angle(self, angle):
        self._angles.append(angle)

    def _build_propers(self):
        periodicities = self.parameters.get_dihedral_periodicities()
        phases = self.parameters.get_dihedral_phases()
        ks = self.parameters.get_dihedral_ks()
        idivfs = self.parameters.get_dihedral_idivfs()

        # idivf is a optional parameter in OpenForceField
        if len(idivfs) == 0:
            for period_by_index in periodicities:
                idivfs.append(dict(zip(period_by_index.keys(),
                                       [1, ] * len(period_by_index.keys()))))

        assert len(periodicities) == len(phases) and \
            len(periodicities) == len(ks) and \
            len(periodicities) == len(idivfs), 'Unconsistent set of ' \
            'OpenForceField\'s torsional parameters. They all should have ' \
            'equal lengths'

        for period_by_index, phase_by_index, k_by_index, idivf_by_index in \
                zip(periodicities, phases, ks, idivfs):

            assert period_by_index.keys() == phase_by_index.keys() and \
                period_by_index.keys() == k_by_index.keys() and \
                period_by_index.keys() == idivf_by_index.keys(), 'Unconsistent ' \
                'torsional parameter indexes. Keys should match.'

            for index in period_by_index.keys():
                atom1_idx, atom2_idx, atom3_idx, atom4_idx = index

                period = period_by_index[index]
                phase = phase_by_index[index]
                k = k_by_index[index]
                idivf = idivf_by_index[index]

                if period and phase and k and idivf:
                    off_proper = OFFProper(atom1_idx=atom1_idx,
                                           atom2_idx=atom2_idx,
                                           atom3_idx=atom3_idx,
                                           atom4_idx=atom4_idx,
                                           periodicity=period,
                                           phase=phase,
                                           k=k,
                                           idivf=idivf)

                    PELE_proper = off_proper.to_PELE()
                    self._add_proper(PELE_proper)
                    self._add_OFF_proper(off_proper)

    def _add_proper(self, proper):
        self._propers.append(proper)

    def _add_OFF_proper(self, proper):
        self._OFF_propers.append(proper)

    def _build_impropers(self):
        periodicities = self.parameters.get_improper_periodicities()
        phases = self.parameters.get_improper_phases()
        ks = self.parameters.get_improper_ks()
        idivfs = self.parameters.get_improper_idivfs()

        # idivf is a optional parameter in OpenForceField
        if len(idivfs) == 0:
            for period_by_index in periodicities:
                idivfs.append(dict(zip(period_by_index.keys(),
                                       [1, ] * len(period_by_index.keys()))))

        assert len(periodicities) == len(phases) and \
            len(periodicities) == len(ks) and \
            len(periodicities) == len(idivfs), 'Unconsistent set of ' \
            'OpenForceField\'s improper parameters. They all should have ' \
            'equal lengths'

        for period_by_index, phase_by_index, k_by_index, idivf_by_index in \
                zip(periodicities, phases, ks, idivfs):

            assert period_by_index.keys() == phase_by_index.keys() and \
                period_by_index.keys() == k_by_index.keys() and \
                period_by_index.keys() == idivf_by_index.keys(), 'Unconsistent ' \
                'torsional parameter indexes. Keys should match.'

            for index in period_by_index.keys():
                atom1_idx, atom2_idx, atom3_idx, atom4_idx = index

                period = period_by_index[index]
                phase = phase_by_index[index]
                k = k_by_index[index]
                idivf = idivf_by_index[index]

                if period and phase and k and idivf:
                    off_improper = OFFImproper(atom1_idx=atom1_idx,
                                               atom2_idx=atom2_idx,
                                               atom3_idx=atom3_idx,
                                               atom4_idx=atom4_idx,
                                               periodicity=period,
                                               phase=phase,
                                               k=k,
                                               idivf=idivf)

                    PELE_improper = off_improper.to_PELE()
                    self._add_improper(PELE_improper)
                    self._add_OFF_improper(off_improper)

    def _add_improper(self, improper):
        self._impropers.append(improper)

    def _add_OFF_improper(self, improper):
        self._OFF_impropers.append(improper)

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
    def forcefield(self):
        return self._forcefield

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def propers(self):
        return self._propers

    @property
    def impropers(self):
        return self._impropers
