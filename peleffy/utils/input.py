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

    def extract_molecule_from_chain(self, chain, rotamer_resolution,
                                    exclude_terminal_rotamers,
                                    allow_undefined_stereo):
        """
        It extracts a peleffy.topology.Molecule object selected by the chain.

        Parameters
        ----------
        chain_id : str
            Chain ID.
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        allow_undefined_stereo : bool
            Whether to allow a molecule with undefined stereochemistry
            to be defined or try to assign the stereochemistry and
            raise a complaint if not possible. Default is False

        Returns
        -------
        molecule : a peleffy.topology.Molecule object
            Selected molecule.
        """
        from peleffy.topology.molecule import Molecule

        # Select which atoms compose this molecule
        atom_ids = [line[6:11].strip() for line in self.pdb_content
                    if line.startswith('HETATM') and line[21:22] == chain]

        # Extract the PDB block of the molecule
        pdb_block = [line for line in self.pdb_content
                     if (line.startswith('HETATM') or line.startswith('CONECT'))
                     and any(a in line for a in atom_ids)]
        return Molecule(pdb_block=''.join(pdb_block),
                        rotamer_resolution=rotamer_resolution,
                        exclude_terminal_rotamers=exclude_terminal_rotamers,
                        allow_undefined_stereo=allow_undefined_stereo)

    def get_hetero_molecules(self, rotamer_resolution=30,
                             exclude_terminal_rotamers=True,
                             allow_undefined_stereo=False):
        """
        It returns a list of peleffy.topology.Molecule objects with all the
        hetero molecules contained in the PDB.

        Returns
        -------
        molecules : list[peleffy.topology.Molecule]
            List of the multiple molecules in the PDB file
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        allow_undefined_stereo : bool
            Whether to allow a molecule with undefined stereochemistry
            to be defined or try to assign the stereochemistry and
            raise a complaint if not possible. Default is False
        """
        chain_ids = set([line[21:22] for line in self.pdb_content
                         if line.startswith('HETATM')])
        molecules = [self.extract_molecule_from_chain(
            chain=chain_id,
            rotamer_resolution=rotamer_resolution,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            allow_undefined_stereo=allow_undefined_stereo)
            for chain_id in chain_ids]

        return molecules

    def get_molecule_from_chain(self, selected_chain, rotamer_resolution=30,
                                exclude_terminal_rotamers=True,
                                allow_undefined_stereo=False):
        """
        It selects a molecule from a chain. It handles the possibles error when
        selecting the chain for a PDB, and if any it returns the molecule as a
        peleffy.topology.Molecule object.

        Parameters
        ----------
        selected_chain : str
            Chain Id.
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        allow_undefined_stereo : bool
            Whether to allow a molecule with undefined stereochemistry
            to be defined or try to assign the stereochemistry and
            raise a complaint if not possible. Default is False

        Returns
        -------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object corresponding to the selected chain.
        """

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
        return self.extract_molecule_from_chain(chain=selected_chain,
                            rotamer_resolution=rotamer_resolution,
                            exclude_terminal_rotamers=exclude_terminal_rotamers,
                            allow_undefined_stereo=allow_undefined_stereo)
