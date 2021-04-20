"""
This module contains a set of classes and methods designed to handle
input data.
"""
from peleffy.utils import Logger


class PDBFile(object):
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

        >>> from peleffy.utils.input import PDBFile

        >>> PDBreader = PDBFile('/path/to/pdb.pdb')
        >>> molecules = PDBreader.get_hetero_molecules()

        Load from a PDB the hetero atom in the L chain.

        >>> from peleffy.utils.input import PDBFile

        >>> PDBreader = PDBFile('/path/to/pdb.pdb')
        >>> molecule  = PDBreader.get_molecule_from_chain(selected_chain = 'L')

        """
        with open(path, 'r') as pdb_file:
            self.pdb_content = pdb_file.readlines()

    def _extract_molecules_from_chain(self, chain, rotamer_resolution,
                                     exclude_terminal_rotamers,
                                     allow_undefined_stereo,
                                     core_constraints):
        """
        It extracts all hetero molecules found in the selected the chain
        of a PDB file.

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
        core_constraints : list[int or str]
            It defines the list of atoms to constrain in the core, thus,
            the core will be forced to contain them. Atoms can be specified
            through integers that match the atom index or strings that
            match with the atom PDB name

        Returns
        -------
        molecules : list[peleffy.topology.Molecule object]
            Selected molecules
        """
        from peleffy.topology.molecule import Molecule

        # Check if there is more than one hetero molecule in the same chain
        residues_ids = set([line[22:26].strip() for line in self.pdb_content
                            if line.startswith('HETATM')
                            and line[21:22] == chain
                            and not line[17:20].strip() == 'HOH'])

        molecules = []
        for residue_id in residues_ids:
            res_name = set([line[17:20].strip() for line in self.pdb_content
                            if line.startswith('HETATM')
                            and line[21:22] == chain
                            and line[22:26].strip() == residue_id])

            # Select which atoms compose this hetero molecule
            atom_ids = [line[6:11].strip() for line in self.pdb_content
                        if line.startswith('HETATM')
                        and line[21:22] == chain
                        and line[22:26].strip() == residue_id]

            # Extract the PDB block of the molecule
            pdb_block = [line for line in self.pdb_content
                         if (line.startswith('HETATM') or
                             line.startswith('CONECT'))
                         and any(' {} '.format(a) in line for a in atom_ids)]

            try:
                molecules.append(
                    Molecule(pdb_block=''.join(pdb_block),
                             rotamer_resolution=rotamer_resolution,
                             exclude_terminal_rotamers=exclude_terminal_rotamers,
                             allow_undefined_stereo=allow_undefined_stereo,
                             core_constraints=core_constraints))
            except Exception as e:
                log = Logger()
                log.warning(' - Skipping {} '.format(list(res_name)[0])
                            + 'from chain {}'.format(chain))
                log.warning('  - The following exception was raised: '
                            + '{}'.format(e))

        return molecules

    def get_hetero_molecules(self, rotamer_resolution=30,
                             exclude_terminal_rotamers=True,
                             allow_undefined_stereo=False,
                             ligand_core_constraints=[],
                             ligand_resname=None):
        """
        It returns a list of peleffy.topology.Molecule objects with all the
        hetero molecules contained in the PDB.

        Parameters
        ----------
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
        ligand_core_constraints : list[int or str]
            It defines the list of atoms to constrain in the core of the
            ligand, thus, the core will be forced to contain them.
            Atoms can be specified through integers that match the atom index
            or strings that match with the atom PDB name
        ligand_resname : str
            Residue name of the ligand. Default is None

        Returns
        -------
        molecules : list[peleffy.topology.Molecule]
            List of the multiple molecules in the PDB file
        """
        chain_ids = set([line[21:22] for line in self.pdb_content
                         if line.startswith('HETATM')
                         and not line[17:20].strip() == 'HOH'])

        # Assign core constraints to the ligand, if specified
        d = {}
        for chain in chain_ids:
            resnames = \
                list(set([line[17:20].strip() for line in self.pdb_content
                          if line.startswith('HETATM') and line[21:22] == chain
                          and not line[17:20].strip() == 'HOH']))
            core_constraints = \
                ligand_core_constraints if ligand_resname in resnames else []
            d[chain] = {'Residues names': resnames,
                        'Core constraints': core_constraints}

        molecules = [self._extract_molecules_from_chain(
            chain=chain_id,
            rotamer_resolution=rotamer_resolution,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            allow_undefined_stereo=allow_undefined_stereo,
            core_constraints=d[chain_id]['Core constraints'])
            for chain_id in chain_ids]

        return sum(molecules, [])

    def get_molecules_from_chain(self, selected_chain, rotamer_resolution=30,
                                 exclude_terminal_rotamers=True,
                                 allow_undefined_stereo=False,
                                 core_constraints=[]):
        """
        It returns all hetero molecule defined in a specific chain from the
        PDB file. It handles the possible errors when selecting the chain
        for a PDB.

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
        molecule : list[peleffy.topology.Molecule]
            The list of peleffy's Molecule object corresponding to the
            selected chain
        """

        chain_ids = set([line[21:22] for line in self.pdb_content
                         if line.startswith('HETATM')
                         and not line[17:20].strip() == 'HOH'])

        all_chain_ids = set([line[21:22] for line in self.pdb_content
                             if line.startswith('ATOM')
                             or line.startswith('HETATM')])

        if not selected_chain in all_chain_ids:
            raise ValueError('The selected chain {} '.format(selected_chain)
                             + 'is not a valid chain for this PDB. Available '
                             + 'chains to select are: {}'.format(chain_ids))

        if not selected_chain in chain_ids and selected_chain in all_chain_ids:
            raise ValueError('The selected chain {} '.format(selected_chain)
                             + 'is not an hetero molecule. Peleffy '
                             + 'is only compatible with hetero atoms.')

        molecules = self._extract_molecules_from_chain(
            chain=selected_chain,
            rotamer_resolution=rotamer_resolution,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            allow_undefined_stereo=allow_undefined_stereo,
            core_constraints=core_constraints)

        return molecules
