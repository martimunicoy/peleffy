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
                                      core_constraints,
                                      only_hetero=True):
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
        only_hetero : bool
            Whether to extract only hetero molecules or not.
            Default is True

        Returns
        -------
        molecules : list[peleffy.topology.Molecule object]
            Selected molecules
        """
        from peleffy.topology.molecule import Molecule

        # Check if there is more than one hetero molecule in the same chain
        residues_ids = set([line[22:26].strip() for line in self.pdb_content
                            if (line.startswith('HETATM') or line.startswith('ATOM') and only_hetero)
                            and line[21:22] == chain
                            and not line[17:20].strip() == 'HOH'])

        molecules = []
        for residue_id in residues_ids:
            res_name = set([line[17:20].strip() for line in self.pdb_content
                            if (line.startswith('HETATM') or line.startswith('ATOM') and only_hetero)
                            and line[21:22] == chain
                            and line[22:26].strip() == residue_id])

            # Select which atoms compose this hetero molecule
            atom_ids = [line[6:11].strip() for line in self.pdb_content
                        if (line.startswith('HETATM') or line.startswith('ATOM') and only_hetero)
                        and line[21:22] == chain
                        and line[22:26].strip() == residue_id]

            # Extract the PDB block of the molecule
            pdb_block = []
            for line in self.pdb_content:

                if ((line.startswith('HETATM') or line.startswith('ATOM') and only_hetero)
                        and line[6:11].strip() in atom_ids):
                    pdb_block.append(line)

                if line.startswith('CONECT'):
                    stripped_line = line.replace("CONECT", "")
                    ids_in_line = [stripped_line[i:i + 5] for i in
                                   range(0, len(stripped_line), 5)]

                    # Strip out whitespaces from ids
                    stripped_ids_in_line = [element.strip() for
                                            element in ids_in_line]

                    if any([atom_id in stripped_ids_in_line for
                            atom_id in atom_ids]):
                        pdb_block.append(line)

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
            Chain Id
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
                         if line.startswith('ATOM')
                         or line.startswith('HETATM')])

        if not selected_chain in chain_ids:
            raise ValueError('The selected chain {} '.format(selected_chain)
                             + 'is not a valid chain for this PDB. Available '
                             + 'chains to select are: {}'.format(chain_ids))

        molecules = self._extract_molecules_from_chain(
            chain=selected_chain,
            rotamer_resolution=rotamer_resolution,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            allow_undefined_stereo=allow_undefined_stereo,
            core_constraints=core_constraints)

        return molecules
    
    def get_molecule_from_residue_id(self, molecules_in_chain, selected_residue_id,
                                     rotamer_resolution=30,
                                     exclude_terminal_rotamers=True,
                                     allow_undefined_stereo=False,
                                     core_constraints=[]):
        """
        It returns all hetero molecule defined in a specific chain from the
        PDB file. It handles the possible errors when selecting the chain
        for a PDB.

        Parameters
        ----------
        molecules_in_chain : list[peleffy.topology.Molecule]
            List of peleffy's Molecule objects to search for the selected
            residue id
        selected_residue_id : int
            Residue Id
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
        molecule : peleffy.topology.Molecule
            The peleffy's Molecule object corresponding to the
            selected residue id
        """
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        # Initialize RDKit toolkit
        rdkit = RDKitToolkitWrapper()

        # Get selected chain
        selected_chain = rdkit.get_chain_id(molecules_in_chain[0])

        if selected_chain is None:
            raise ValueError('The selected chain is undefined')

        residue_ids = set([line[22:26].strip() for line in self.pdb_content
                           if line.startswith('ATOM') or line.startswith('HETATM')
                           and line[21:22] == selected_chain
                           and not line[17:20].strip() == 'HOH'])
        
        if len(residue_ids) == 0:
            raise ValueError('The selected chain {} '.format(selected_chain)
                             + 'does not contain any residues')
                
        for molecule, residue_id in zip(molecules_in_chain, residue_ids):
            if residue_id == str(selected_residue_id):
                return molecule
        
        raise ValueError('The selected chain {} '.format(selected_chain)
                         + 'does not contain a residue with id '
                         + '{}'.format(selected_residue_id))
    
    def get_covalent_molecules_to_residue_id(self, molecule, target_molecules):
        """
        It returns molecules defined in the input list that are
        covalently bonded to the selected molecule.

        Parameters
        ----------
        molecule : peleffy.topology.Molecule
            The peleffy's Molecule object corresponding to the selected
            residue id
        target_molecules : list[peleffy.topology.Molecule]
            List of peleffy's Molecule objects to search for the selected
            residue id
        
        Returns
        -------
        covalent_molecules : list[peleffy.topology.Molecule]
            The list of peleffy's Molecule objects corresponding to the
            selected residue id
        """
        import numpy as np
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        # Initialize RDKit toolkit
        rdkit = RDKitToolkitWrapper()

        # Extract coordinates of the selected residue
        molecule_coordinates = rdkit.get_coordinates(molecule)
                
        # Iterate over all molecules to find the covalent ones
        covalent_molecules = []
        for target_molecule in target_molecules:
            if rdkit.get_residue_id(target_molecule) == rdkit.get_residue_id(molecule):
                continue
            
            # Extract coordinates of the molecule
            target_molecule_coordinates = rdkit.get_coordinates(target_molecule)

            # Check if the molecule is covalently bonded to the selected residue
            distances = np.linalg.norm(molecule_coordinates[:, np.newaxis, :] - target_molecule_coordinates[np.newaxis, :, :], axis=2)
            lowest_distance = np.min(distances)

            if lowest_distance <= 1.5:
                covalent_molecules.append(target_molecule)
        
        return covalent_molecules

    @property
    def is_complex(self):
        """
        Check whether the PDB fetched corresponds to a protein-ligand complex 
        or not.

        Returns
        -------
        is_complex : bool
            True if it is a protein-ligand complex.
        """

        if any(line.startswith('ATOM') for line in self.pdb_content) and \
           any(line.startswith('HETATM') for line in self.pdb_content):
            return True
        else:
            return False

    @staticmethod
    def is_unique(molecules):
        """
        Check whether a list of molecules contains only one or multiple 
        elements. 

        Parameters
        ----------
        molecule : list[peleffy.topology.Molecule]
            A list of peleffy's Molecule object.

        Returns
        -------
        is_unique : bool
            True if it only contains one molecule.
        """
        return len(molecules) == 1
