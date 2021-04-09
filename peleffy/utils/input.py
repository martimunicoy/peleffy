"""
This module contains a set of classes and methods designed to handle
input data.
"""


class PDB(object):
    """
    It handles an input PDB file and allows the extraction of multiple molecules
    as peleffy.topology.Molecule objects.
    """

    def __init__(self, path):
        """
        It initializes a PDB object through a PDB file.

        Parameters
        ----------
        path : str
            The path to the PDB with the molecules structures.

        Examples
        --------

        Load all the hetero molecules from a PDB

        >>> from peleffy.utils.input import PDB

        >>> PDBreader = PDB('/path/to/pdb.pdb')
        >>> molecules = PDBreader.get_hetero_molecules()

        Load from a PDB the hetero atom in the L chain.

        >>> from peleffy.utils.input import PDB

        >>> PDBreader = PDB('/path/to/pdb.pdb')
        >>> molecule  = PDBreader.get_molecule_from_chain(selected_chain = 'L')

        """
        self.pdb_content = open(path, 'r').readlines()

    def extract_molecule_from_chain(self, chain_id):
        """
        It extracts a peleffy.topology.Molecule object selected by the chain.

        Parameters
        ----------
        chain_id : str
            Chain ID.

        Returns
        -------
        molecule : a peleffy.topology.Molecule object
            Selected molecule.
        """
        from peleffy.topology.molecule import Molecule

        # Select which atoms compose this molecule
        atom_ids = [line[6:11].strip() for line in self.pdb_content
                    if line.startswith('HETATM') and line[21:22] == chain_id]

        # Extract the PDB block of the molecule
        pdb_block = [line for line in self.pdb_content
                     if (line.startswith('HETATM') or line.startswith('CONECT'))
                     and any(a in line for a in atom_ids)]

        return Molecule(pdb_block=''.join(pdb_block))

    def get_hetero_molecules(self):
        """
        It returns a list of peleffy.topology.Molecule objects with all the
        hetero molecules contained in the PDB.

        Returns
        -------
        molecules : list[peleffy.topology.Molecule]
            List of the multiple molecules in the PDB file
        """
        chain_ids = set([line[21:22] for line in self.pdb_content
                         if line.startswith('HETATM')])
        molecules = [self.get_molecule_from_chain(selected_chain=chain_id) for
                     chain_id in chain_ids]

        return molecules

    def get_molecule_from_chain(self, selected_chain):
        """
        It selects a molecule from a chain. It handles the possibles error when
        selecting the chain for a PDB, and if any it returns the molecule as a
        peleffy.topology.Molecule object.

        Parameters
        ----------
        selected_chain : str
            Chain Id.

        Returns
        -------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object corresponding to the selected chain.
        """
        if selected_chain is None:
            raise ValueError('The input PDB has multiple molecules. A chain' +
                             'has to be selected.')

        chain_ids = set([line[21:22] for line in self.pdb_content
                         if line.startswith('HETATM')])
        all_chain_ids = set([line[21:22] for line in self.pdb_content
                             if line.startswith('ATOM')
                             or line.startswith('HETATM')])
        if not selected_chain in all_chain_ids:
            raise ValueError('The selected chain {}'.format(selected_chain) +
                             ' is not a valid chain for this PDB. Available' +
                             ' chains to select are: {}'.format(chain_ids))
        if not selected_chain in chain_ids and selected_chain in all_chain_ids:
            raise ValueError('The selected chain {}'.format(selected_chain) +
                             ' is not a hetero molecule. Peleffy' +
                             ' is only compatible with hetero atoms.')
        return self.extract_molecule_from_chain(chain_id=selected_chain)
