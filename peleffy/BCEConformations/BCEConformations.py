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
from peleffy.template import impact
from peleffy.topology import molecule
from peleffy.utils.toolkits import RDKitToolkitWrapper

class BCEConformations(object):
    """
    A class to produce a library of conformations from the output
    of the BCE server
    """
    def __init__(self, topology_obj, bce_output_path, verbose=False):
        """
        Initializes a BCEConformations object

        Parameters
        ----------
        topology_obj : An peleffy.topology.Topology
            A Topology object that contains the ligand's information
        bce_output_path: str
            Path where the output from the BCE server is stored

        """
        self._topology = topology_obj
        self._molecule = topology_obj.molecule
        self.bce_path = bce_output_path
        self.conformations_library = {}
        self.verbose = verbose

    def calculate(self):
        """
        Calculate conformation library from the bce output
        """
        logger = Logger()
        logger.info(' - Calculating conformation library')
        if not self.verbose:
            logger.set_level('WARNING')
        self._calculate_all_conformations()

    def _calculate_all_conformations(self):
        clusters = sorted(glob.glob(os.path.join(self.bce_path, "CLUSTERS", "CL*", "cluster*.min.imaged.pdb")))
        if not clusters:
            raise ValueError("Path to the BCE output does not contain a CLUSTERS folder, please check if the path is correct!")
        clusters_order = self.order_clusters_min_distances(clusters)
        ordered_clusters = [clusters[x] for x in clusters_order]
        for cluster in ordered_clusters:
            self.calculate_cluster_offsets(cluster)


    def calculate_cluster_offsets(self, cluster_pdb):
        """
        Calculate dihedral angles from pdb

        Parameters
        ----------
        cluster_pdb: str
            Path to the cluster representative conformation
        """
        rdkit_wrapper = RDKitToolkitWrapper()
        conformation_offsets = []
        # use the input molecule as template since the cluster structures
        # probably will not have proper stereochemistry
        mol = molecule.Molecule(cluster_pdb, connectivity_template=self._molecule.rdkit_molecule)
        cluster_coordinates = mol.get_conformer()
        topology_to_cluster = {i: x for i, x in enumerate(rdkit_wrapper.get_substruct_match(mol, self._molecule))}
        cluster_to_topology = {v: k for k, v in topology_to_cluster.items()}
        template = impact.Impact(self._topology)
        root_coordinates = cluster_coordinates[topology_to_cluster[find_index_root(template.topology, self._topology)]]
        for i, coords in enumerate(cluster_coordinates):
            conformation_offsets.append([self._topology.atoms[cluster_to_topology[i]].PDB_name]+(coords-root_coordinates).tolist())
        self.conformations_library[cluster_pdb] = conformation_offsets

    def save(self, output_path):
        """
        Save the conformation library

        Parameters
        ----------
        output_path: str
            Where to save the conformation library
        """
        with open(output_path, "w") as fw:
            fw.write("* CONFORMATIONS LIBRARY FILE\n")
            fw.write('* File generated with peleffy-{}\n'.format( peleffy.__version__))
            for conformation_file, conformation_values in self.conformations_library.items():
                fw.write("* File: {:s}\n".format(conformation_file))
                fw.write("{:s} {:d} {:d}\n".format(self._molecule.tag.upper(), len(conformation_values), len(self.conformations_library)))
                for values in conformation_values:
                    fw.write("{:s} {:5f} {:5f} {:5f}\n".format(*values))
                fw.write("ENDCONFORMATION\n")
            fw.write("END\n")

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


def find_index_root(sorted_topology, topology):
    atom_name = sorted_topology.atoms[0].PDB_name
    for i, atom in enumerate(topology.atoms):
        if atom.PDB_name == atom_name:
            return i
    return None
