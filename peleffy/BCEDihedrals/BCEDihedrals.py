"""
This module contains classes and functions involved in the
creation of a library of dihedral angles from the output of the
BCE server
"""
import os
import glob
import numpy as np
import networkx as nx
import peleffy
from peleffy.utils import Logger
from peleffy.topology import molecule
from peleffy.utils.toolkits import RDKitToolkitWrapper

class BCEDihedrals(object):
    """
    A class to produce a library of dihedral angles from the output
    of the BCE server
    """
    def __init__(self, topology_obj, bce_output_path, mode="all_dihedrals", verbose=False):
        """
        Initializes a BCEDihedrals object

        Parameters
        ----------
        topology_obj : An peleffy.topology.Topology
            A Topology object that contains the ligand's information
        bce_output_path: str
            Path where the output from the BCE server is stored
        mode: str
            Whether to extract all dihedrals or only those marked as flexible

        """
        self._topology = topology_obj
        self._molecule = topology_obj.molecule
        self.bce_path = bce_output_path
        self.dihedral_library = {}
        self.mode = mode
        self.verbose = verbose

    def calculate(self):
        """
        Calculate dihedrals library from the bce output
        """
        logger = Logger()
        logger.info(' - Calculating dihedral library')
        logger.info('   - with mode: {}'.format(self.mode))
        if not self.verbose:
            logger.set_level('WARNING')
        if self.mode == "all_dihedrals":
            self._calculate_all_dihedrals()
        elif self.mode == "flexible_dihedrals":
            self._calculate_flexible_dihedrals()

    def _calculate_all_dihedrals(self):
        clusters = sorted(glob.glob(os.path.join(self.bce_path, "CLUSTERS", "CL*", "cluster*.min.imaged.pdb")))
        if not clusters:
            raise ValueError("Path to the BCE output does not contain a CLUSTERS folder, please check if the path is correct!")
        clusters_order = self.order_clusters_min_distances(clusters)
        dihedrals = self.list_all_dihedrals()
        ordered_clusters = [clusters[x] for x in clusters_order]
        for cluster in ordered_clusters:
            self.calculate_cluster_angles(cluster, dihedrals, match_indexes=False)

    def _calculate_flexible_dihedrals(self):
        # here we rely on the output of the bce server to not change its
        # structrue in the future
        dihedral_file = glob.glob(os.path.join(self.bce_path, "DIHEDRALS", "*_new.ndx"))[0]
        dihedrals = read_dihedral_info_file(dihedral_file)
        clusters = sorted(glob.glob(os.path.join(self.bce_path, "CLUSTERS", "CL*", "cluster*.min.imaged.pdb")))
        for cluster in clusters:
            self.calculate_cluster_angles(cluster, dihedrals)


    def calculate_cluster_angles(self, cluster_pdb, dihedral_list, match_indexes=True):
        """
        Calculate dihedral angles from pdb

        Parameters
        ----------
        cluster_pdb: str
            Path to the cluster representative conformation
        dihedral_list: list
            List of the tuples containing the atoms that form the dihedrals
        match_indexes: bool
            Whether to use the atom indices from the dihedral list or match to
            the cluster structure before
        """
        rdkit_wrapper = RDKitToolkitWrapper()
        pdb_dihedrals = []
        # use the input molecule as template since the cluster structures
        # probably will not have proper stereochemistry
        mol = molecule.Molecule(cluster_pdb, connectivity_template=self._molecule.rdkit_molecule)
        # we use substructure matching to ensure that the indices in the
        # clusters pdb and the input ligand are the same
        if match_indexes:
            rename_dict = {i: x for i, x in enumerate(rdkit_wrapper.get_substruct_match(self._molecule, mol))}
        else:
            rename_dict = {i: x for i, x in enumerate(rdkit_wrapper.get_substruct_match(mol, self._molecule))}
        for dihedral in dihedral_list:
            if match_indexes:
                names = [self._topology.atoms[rename_dict[atom]].PDB_name for atom in dihedral]
                angle = rdkit_wrapper.get_dihedral(mol, *dihedral, units="radians")
            else:
                names = [self._topology.atoms[atom].PDB_name for atom in dihedral]
                dihedral_matched = tuple([rename_dict[index] for index in dihedral])
                angle = rdkit_wrapper.get_dihedral(mol, *dihedral_matched, units="radians")

            pdb_dihedrals.append(names+[angle])
        self.dihedral_library[cluster_pdb] = pdb_dihedrals

    def save(self, output_path):
        """
        Save the dihedral library

        Parameters
        ----------
        output_path: str
            Where to save the dihedral library
        """
        with open(output_path, "w") as fw:
            fw.write("* DIHEDRAL LIBRARY FILE\n")
            fw.write('* File generated with peleffy-{}\n'.format( peleffy.__version__))
            for dihedral_file, dihedral_values in self.dihedral_library.items():
                fw.write("* File: {:s}\n".format(dihedral_file))
                fw.write("{:s} {:d} {:d}\n".format(self._molecule.tag.upper(), len(dihedral_values), len(self.dihedral_library)))
                for dihedral in dihedral_values:
                    fw.write("{:s} {:s} {:s} {:s} {:5f}\n".format(*dihedral))
                fw.write("ENDDIHEDRALS\n")
            fw.write("END")

    def list_all_dihedrals(self):
        dihedrals = []
        seen_dihedrals = set()
        for proper in self._topology.propers:
            dihedral_indexes = (proper.atom1_idx, proper.atom2_idx, proper.atom3_idx, proper.atom4_idx)
            if dihedral_indexes in seen_dihedrals:
                continue
            seen_dihedrals.add(dihedral_indexes)
            dihedrals.append(list(dihedral_indexes))
        return dihedrals


    def order_clusters_min_distances(self, clusters):
        """
        Order a list of structures in a way that jumps to the next structure
        involves the minimum change

        Parameters
        ----------

        clusters: list
            Path to the structures
        """
        rdkit_wrapper = RDKitToolkitWrapper()
        distances = np.zeros((len(clusters), len(clusters)))
        cluster_molecules = []
        for cluster in clusters:
            # use the input molecule as template since the cluster structures
            # probably will not have proper stereochemistry
            cluster_molecules.append(molecule.Molecule(cluster, connectivity_template=self._molecule.rdkit_molecule))
        for i, cluster_mol in enumerate(cluster_molecules):
            for j, cluster_mol_2 in enumerate(cluster_molecules[i+1:], start=i+1):
                rmsd_value = rdkit_wrapper.get_rmsd(cluster_mol, cluster_mol_2)
                distances[i, j] = rmsd_value
                distances[j, i] = rmsd_value
        return find_optimal_path_from_matrix(distances)


def find_optimal_path_from_matrix(distances):
    """
    Find a Minimum Spanning Tree from a distance matrix. It uses networkx to
    create the graph and the MST

    Parameters
    ----------

    Returns
    -------
    min_path : list
        Nodes on the graph that minimise heuristically the total distance
    """
    graph = nx.from_numpy_matrix(distances)
    min_dist = 1e6
    min_path = None
    for node in graph:
        dist, path = find_heuristic_path(graph.copy(), distances, node)
        if dist < min_dist:
            min_dist = dist
            min_path = path
    return min_path

def find_heuristic_path(graph, distances, start_node):
    n = distances.shape[0]
    path = [start_node]
    dist = 0
    while len(path) < n:
        node = path[-1]
        new_neighbor = min(graph[node].items(), key=lambda x: x[1]["weight"])
        graph.remove_node(node)
        path.append(new_neighbor[0])
        dist += new_neighbor[1]['weight']
    dist += distances[path[-1], start_node]
    return dist, path

def read_dihedral_info_file(file_name):
    """
    Read an ndx file containing the classification of all dihedrals and return
    the flexible ones

    Parameters
    ----------
    file_name: str
        Path of the dihedrals file

    Returns
    -------
    dihedrals : list
        List of atom indices correponding to the flexible dihedrals
    """
    read_dihedral = False
    dihedrals = []
    with open(file_name) as f:
        for line in f:
            # reached the block of all dihedrals, we can stop reading
            if line.startswith("[ ALL ]"):
                return dihedrals
            if line.startswith("["):
                name = line.strip("\n[] ")
                read_dihedral = name.startswith("F")
            else:
                if read_dihedral:
                    # the indices of the bce server are 1-based, while our
                    # indices are 0-based
                    dihedrals.append([int(x)-1 for x in line.rstrip().split()])
                    read_dihedral = False
    return dihedrals
