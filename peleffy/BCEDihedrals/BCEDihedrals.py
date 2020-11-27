"""
This module contains classes and functions involved in the
creation of a library of dihedral angles from the output of the
BCE server
"""
import os
import glob
import peleffy
from peleffy.utils import Logger
from peleffy.topology import molecule
from peleffy.utils.toolkits import RDKitToolkitWrapper

class BCEDihedrals(object):
    """
    A class to produce a library of dihedral angles from the output
    of the BCE server
    """
    def __init__(self, topology_obj, bce_output_path):
        """
        Initializes a BCEDihedrals object

        Parameters
        ----------
        molecule_obj : An peleffy.topology.Topology
            A Topology object that contains the ligand's information

        bce_output_path: str
            Path where the output from the BCE server is stored
        """
        self._topology = topology_obj
        self._molecule = topology_obj.molecule
        self.bce_path = bce_output_path
        self.dihedral_library = {}

    def calculate(self):
        """
        Calculate dihedrals library from the bce output
        """
        logger = Logger()
        logger.info(' - Calculating dihedral library')
        # here we rely on the output of the bce server to not change its
        # structrue in the future
        dihedral_file = glob.glob(os.path.join(self.bce_path, "DIHEDRALS", "*_new.ndx"))[0]
        dihedrals = read_dihedral_info_file(dihedral_file)
        clusters = glob.glob(os.path.join(self.bce_path, "CLUSTERS", "CL*", "cluster*.min.imaged.pdb"))
        for cluster in clusters:
            self.calculate_cluster_angles(cluster, dihedrals)


    def calculate_cluster_angles(self, cluster_pdb, dihedral_list):
        """
        Calculate dihedral angles from pdb

        Parameters
        ----------
        cluster_pdb: str
            Path to the cluster representative conformation

        dihedral_list: list
            List of the tuples containing the atoms that form the dihedrals
        """
        rdkit_wrapper = RDKitToolkitWrapper()
        pdb_dihedrals = []
        # use the input molecule as template since the cluster structures
        # probably will not have proper stereochemistry
        mol = molecule.Molecule(cluster_pdb, connectivity_template=self._molecule.rdkit_molecule)
        # we use substructure matching to ensure that the indices in the
        # clusters pdb and the input ligand are the same
        rename_dict = {i: x for i, x in enumerate(rdkit_wrapper.get_substruct_match(self._molecule, mol))}
        for dihedral in dihedral_list:
            # names = [self._molecule.get_pdb_atom_names()[rename_dict[atom]] for atom in dihedral]
            names = [self._topology.atoms[rename_dict[atom]].PDB_name for atom in dihedral]
            angle = rdkit_wrapper.get_dihedral(mol, *dihedral, units="radians")
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

def read_dihedral_info_file(file_name):
    """
    Read an ndx file containing the classification of all dihedrals and return
    the flexible ones

    Parameters
    ----------
    file_name: str
        Path of the dihedrals file
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
