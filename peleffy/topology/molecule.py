"""
This module handles all classes and functions related with molecular
representations.
"""


from .rotamer import MolecularGraph, MolecularGraphWithConstrainedCore
from peleffy.utils.toolkits import (RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper)
from peleffy.utils import Logger


class Molecule(object):
    """
    It represent wraps up all the tools to parameterize a molecule with
    the OpenForceField toolkit for PELE.
    """

    def __init__(self, path=None, smiles=None, pdb_block=None,
                 rotamer_resolution=30,exclude_terminal_rotamers=True, name='',
                 tag='UNK', connectivity_template=None, core_constraints=[],
                 allow_undefined_stereo=False, fix_pdb=True):
        """
        It initializes a Molecule object through a PDB file or a SMILES
        tag.

        Parameters
        ----------
        path : str
            The path to a PDB with the molecule structure
        smiles : str
            The smiles tag
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        name : str
            The molecule name
        tag : str
            The molecule tag. It must be a 3-character string
        connectivity_template : an rdkit.Chem.rdchem.Mol object
            A molecule represented with RDKit to use when assigning the
            connectivity of this Molecule object
        core_constraints : list[int or str]
            It defines the list of atoms to constrain in the core, thus,
            the core will be forced to contain them. Atoms can be specified
            through integers that match the atom index or strings that
            match with the atom PDB name
        allow_undefined_stereo : bool
            Whether to allow a molecule with undefined stereochemistry
            to be defined or try to assign the stereochemistry and
            raise a complaint if not possible. Default is False
        fix_pdb : bool
            Activates or deactivate the PDB fixer that is executed
            prior parsing it

        Examples
        --------

        Load a molecule from a PDB file and parameterize it with Open
        Force Field

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')
        >>> molecule.parameterize('openff_unconstrained-1.2.0.offxml')

        Load a molecule using a SMILES tag and parameterize it with Open
        Force Field

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule(smiles='Cc1ccccc1')

        Load a molecule usign a PDB file (without connectivity) and assign
        the missing connectivity from an RDKit template (e.g. obtained
        from qcportal and the Open Force Field Toolkit)

        >>> import qcportal as ptl
        >>> from openff.toolkit.topology import Molecule as OFFMolecule

        >>> ds = client.get_collection('OptimizationDataset',
                                       'Kinase Inhibitors: WBO Distributions')
        >>> entry = ds.get_entry(ds.df.index[0])
        >>> mol_record = OFFMolecule.from_qcschema(entry)
        >>> template = mol_record.to_rdkit()

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('PDB_without_connectivity.pdb',
                                template=template)

        Load a molecule and create its rotamer library template with
        a core constraint

        >>> from peleffy.topology import Molecule
        >>> from peleffy.topology import RotamerLibrary

        >>> molecule = Molecule(smiles='CCCC', name='butane', tag='BUT',
                                exclude_terminal_rotamers=False,
                                core_constraints=[0, ])

        >>> rotamer_library = RotamerLibrary(mol)
        >>> rotamer_library.to_file('butz')

        Deactivate OpenFF Toolkit's stereochemistry checking when
        loading a molecule

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule_with_undefined_stereochemistry.pdb',
                                allow_undefined_stereo=True)

        Display the molecular representation in a Jupyter Notebook

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule(smiles='CCCC', name='butane', tag='BUT')

        >>> display(molecule)

        """
        self._name = name
        self._tag = tag
        self._rotamer_resolution = rotamer_resolution
        self._exclude_terminal_rotamers = exclude_terminal_rotamers
        self._connectivity_template = connectivity_template
        self._core_constraints = core_constraints
        self._allow_undefined_stereo = allow_undefined_stereo
        self._fix_pdb = fix_pdb

        # Deactivate OpenForceField toolkit warnings
        import logging
        logging.getLogger().setLevel(logging.ERROR)
        if isinstance(path, str):
            from pathlib import Path
            extension = Path(path).suffix
            extension = extension.strip('.')
            if extension == 'pdb':
                self._initialize_from_pdb(path)
            else:
                raise ValueError(
                    '{} is not a valid extension'.format(extension))
        elif isinstance(smiles, str):
            self._initialize_from_smiles(smiles)
        elif isinstance(pdb_block, str):
            self._initialize_from_pdb_block(pdb_block)
        else:
            self._initialize()

        self._build_rotamers()

    def _initialize(self):
        """It initializes an empty molecule."""
        self._rdkit_molecule = None
        self._off_molecule = None
        self._rotamers = None
        self._graph = None

    def _pdb_checkup(self, path):
        """
        Safety check for PDB files in order to properly handle exceptions
        related with its format prior running the parser.

        Parameters
        ----------
        path : str
            The path to a PDB with the molecule structure
        """

        # Parse PDB file
        atom_id, res_name, res_id = ([] for i in range(3))
        connectivity = False
        with open(path) as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom_id.append(line[12:16])
                    res_name.append(line[17:20])
                    res_id.append(line[22:26])
                if line.startswith('CONECT'):
                    connectivity = True

        # Handle exceptions related with the PDB file format
        if not res_id[:-1] == res_id[1:]:
            raise Exception(
                'A single ligand with immutable residue ids is expected')
        if not res_name[:-1] == res_name[1:]:
            raise Exception(
                'A single ligand with immutable residue names is expected')
        if not len(atom_id) == len(set(atom_id)):
            raise Exception(
                'Ligand in input PDB has no unique atom names')
        if not connectivity and self.connectivity_template is None:
            log = Logger()
            log.warning(
                "Warning: input PDB has no information about the "
                + "connectivity and this could result in an unexpected "
                + "bond assignment")

    def _read_and_fix_pdb(self, path):
        """
        It reads the input PDB file returns the corresponding PDB block.
        It also applies some modifications, in case it requires some
        fixing prior running the parser.

        Parameters
        ----------
        path : str
            The path to a PDB with the molecule structure

        Returns
        -------
        pdb_block : str
            The corresponding PDB block, with applied fixes if required
        """
        log = Logger()

        # Skip PDB fixing if it has been deactivated
        if not self.fix_pdb:
            with open(path) as pdb_file:
                pdb_block = pdb_file.read()

            return pdb_block

        # Fix PDB
        missing_element = False
        any_fail = False
        pdb_block = ''
        with open(path) as pdb_file:
            for line in pdb_file:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    if len(line) < 78 or line[76:78] == '  ':
                        missing_element = True
                        atom_name = line[12:16]
                        # Try to infer element from atom name
                        inferred_element = ''.join([c for c in atom_name
                                                    if not c.isdigit()
                                                    and c != ' '])

                        # Format properly the element identifier
                        if len(inferred_element) == 1:
                            inferred_element = inferred_element.upper()
                        elif len(inferred_element) == 2:
                            inferred_element = inferred_element[0].upper() + \
                                inferred_element[1].lower()
                        else:
                            # We were expecting an element identifier of 1 or 2 chars
                            any_fail = True
                            break

                        # Remove line breaks, if any
                        line = line.strip()

                        # Fill a short line with white spaces
                        while(len(line) < 79):
                            line += ' '

                        # Add element to line (right-justified)
                        line = line[:76] + '{:>2s}'.format(inferred_element) \
                            + line[79:] + '\n'

                pdb_block += line

        if missing_element:
            log.warning(
                "Warning: input PDB has no information about atom "
                + "elements and they were inferred from atom names. "
                + "Please, verify that the resulting elements are "
                + "correct")

        if any_fail:
            log.error("Error: PDB could not be fixed")
            with open(path) as pdb_file:
                pdb_block = pdb_file.read()

        return pdb_block

    def _initialize_from_pdb(self, path):
        """
        It initializes a molecule with the molecule structure read from
        a PDB file.

        Parameters
        ----------
        path : str
            The path to a PDB with the molecule structure
        """
        logger = Logger()
        logger.info(' - Initializing molecule from PDB')
        self._initialize()

        # Validate PDB
        self._pdb_checkup(path)

        # Read and fix PDB
        pdb_block = self._read_and_fix_pdb(path)

        logger.info('   - Loading molecule from RDKit')
        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_pdb_block(pdb_block)

        # Use RDKit template, if any, to assign the connectivity to
        # the current Molecule object
        if self.connectivity_template is not None:
            logger.info('   - Assigning connectivity from template')
            rdkit_toolkit.assign_connectivity_from_template(self)

        if not self.allow_undefined_stereo:
            # RDKit must generate stereochemistry specifically from 3D coords
            logger.info('   - Assigning stereochemistry from 3D coordinates')
            rdkit_toolkit.assign_stereochemistry_from_3D(self)

        # Set molecule name according to PDB name
        if self.name == '':
            from pathlib import Path
            name = Path(path).stem
            logger.info('   - Setting molecule name to \'{}\''.format(name))
            self.set_name(name)

        # Set molecule tag according to PDB's residue name
        if self.tag == 'UNK':
            tag = rdkit_toolkit.get_residue_name(self)
            logger.info('   - Setting molecule tag to \'{}\''.format(tag))
            self.set_tag(tag)

        logger.info('   - Representing molecule with the Open Force Field '
                    + 'Toolkit')
        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        self._off_molecule = openforcefield_toolkit.from_rdkit(self)

    def _initialize_from_smiles(self, smiles):
        """
        It initializes a molecule from a SMILES tag.

        Parameters
        ----------
        smiles : str
            The SMILES tag to construct the molecule structure with
        """
        logger = Logger()
        logger.info(' - Initializing molecule from a SMILES tag')
        self._initialize()

        logger.info('   - Loading molecule from RDKit')
        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_smiles(smiles)

        # TODO not sure if stereochemistry assignment from 3D is still necessary
        # RDKit must generate stereochemistry specifically from 3D coords
        # rdkit_toolkit.assign_stereochemistry_from_3D(self)

        # Set molecule name according to the SMILES tag
        if self.name == '':
            logger.info('   - Setting molecule name to \'{}\''.format(smiles))
            self.set_name(smiles)

        logger.info('   - Representing molecule with the Open Force Field '
                    + 'Toolkit')
        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        self._off_molecule = openforcefield_toolkit.from_rdkit(self)

    def _initialize_from_pdb_block(self, pdb_block):
        """
        It initializes a molecule with the molecule structure fetch in a PDB
        block.

        Parameters
        ----------
        pdb_block : str
            PDB block with the molecule structure
        """
        logger = Logger()
        logger.info(' - Initializing molecule from PDB')
        self._initialize()

        logger.info('   - Loading molecule from RDKit')
        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_pdb_block(pdb_block)

        # Use RDKit template, if any, to assign the connectivity to
        # the current Molecule object
        if self.connectivity_template is not None:
            logger.info('   - Assigning connectivity from template')
            rdkit_toolkit.assign_connectivity_from_template(self)

        if not self.allow_undefined_stereo:
            # RDKit must generate stereochemistry specifically from 3D coords
            logger.info('   - Assigning stereochemistry from 3D coordinates')
            rdkit_toolkit.assign_stereochemistry_from_3D(self)

        # Set molecule tag according to PDB's residue name
        if self.tag == 'UNK':
            tag = rdkit_toolkit.get_residue_name(self)
            logger.info('   - Setting molecule tag to \'{}\''.format(tag))
            self.set_tag(tag)

        logger.info('   - Representing molecule with the Open Force Field '
                    + 'Toolkit')
        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        self._off_molecule = openforcefield_toolkit.from_rdkit(self)

    def _build_rotamers(self):
        """It builds the rotamers of the molecule."""
        logger = Logger()
        if self.off_molecule and self.rdkit_molecule:
            logger.info(' - Generating rotamer library')

            if len(self.core_constraints) != 0:
                self._graph = MolecularGraphWithConstrainedCore(
                    self, self.core_constraints)
                if len(self.core_constraints) == 1:
                    logger.info('   - Core forced to contain atom: '
                                + self._graph.constraint_names[0])
                else:
                    logger.info('   - Core forced to contain atoms: '
                                + ', '.join(atom_name.strip() for atom_name
                                            in self._graph.constraint_names))
            else:
                self._graph = MolecularGraph(self)
                logger.info('   - Core set to the center of the molecule')

            self._rotamers = self._graph.get_rotamers()

    def set_name(self, name):
        """
        It sets the name of the molecule.

        Parameters
        ----------
        name : str
            The name to set to the molecule
        """
        assert isinstance(name, str), 'Invalid type for a name, it must be ' \
            + 'a string'

        self._name = name

    def set_tag(self, tag):
        """
        It sets the tag of the molecule. It must be a 3-character string.

        Parameters
        ----------
        tag : str
            The tag to set to the molecule. It must be a 3-character string
        """
        # Some previous checks
        assert len(tag) == 3, 'Invalid tag length, it must be a ' \
            + '3-character string'
        assert isinstance(tag, str), 'Invalid type for a tag, it must be ' \
            + 'a string'

        self._tag = tag.upper()

    def get_pdb_atom_names(self):
        """
        It returns the PDB atom names of all the atoms in the molecule.

        Returns
        -------
        pdb_atom_names : list[str]
            The PDB atom names of all the atoms in this Molecule object
        """
        rdkit_toolkit = RDKitToolkitWrapper()

        return rdkit_toolkit.get_atom_names(self)

    def to_pdb_file(self, path):
        """
        It writes the molecule to a PDB file.

        Parameters
        ----------
        path : str
            Path to write to
        """
        rdkit_toolkit = RDKitToolkitWrapper()
        rdkit_toolkit.to_pdb_file(self, path)

    def set_conformer(self, conformer):
        """
        It sets a new conformer to the molecule.

        Parameters
        ----------
        conformer : an RDKit.Chem.rdchem.Conformer object
            The conformer to set to the molecule
        """
        rdkit_toolkit = RDKitToolkitWrapper()
        rdkit_toolkit.set_conformer(self, conformer)

    def get_conformer(self):
        """
        It gets the current conformer of the molecule.

        Returns
        -------
        conformer : numpy.ndarray
            The array of 3D coordinates of all the atoms in the molecule
        """
        rdkit_toolkit = RDKitToolkitWrapper()

        return rdkit_toolkit.get_coordinates(self)

    @staticmethod
    def from_rdkit(rdkit_molecule, rotamer_resolution=30,
                   exclude_terminal_rotamers=True, name='', tag='UNK',
                   connectivity_template=None, core_constraints=[],
                   allow_undefined_stereo=False):
        """
        It initializes and returns a peleffy Molecule representation
        from an RDKit molecular representation.

        Parameters
        ----------
        rdkit_molecule : an rdkit.Chem.rdchem.Mol object
            The RDKit's Molecule object to use to initialize a peleffy
            Molecule object
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        name : str
            The molecule name
        tag : str
            The molecule tag. It must be a 3-character string
        connectivity_template : an rdkit.Chem.rdchem.Mol object
            A molecule represented with RDKit to use when assigning the
            connectivity of this Molecule object
        core_constraints : list[int or str]
            It defines the list of atoms to constrain in the core, thus,
            the core will be forced to contain them. Atoms can be specified
            through integers that match the atom index or strings that
            match with the atom PDB name
        allow_undefined_stereo : bool
            Whether to allow a molecule with undefined stereochemistry
            to be defined or try to assign the stereochemistry and
            raise a complaint if not possible. Default is False

        Returns
        -------
        molecule : an peleffy.topology.Molecule
            The resulting peleffy's Molecule object

        Examples
        --------

        Load a molecule from an RDKit molecular representation

        >>> from rdkit import Chem

        >>> rdkit_molecule = Chem.MolFromPDBFile(pdb_path, removeHs=False)

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule.from_rdkit(rdkit_molecule)

        """
        molecule = Molecule(
            rotamer_resolution=30,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            name=name, tag=tag,
            connectivity_template=connectivity_template,
            core_constraints=core_constraints,
            allow_undefined_stereo=allow_undefined_stereo)

        logger = Logger()

        logger.info(' - Initializing molecule from an RDKit '
                    + 'molecular representation')
        molecule._initialize()
        molecule._rdkit_molecule = rdkit_molecule

        logger.info('   - Representing molecule with the Open Force Field '
                    + 'Toolkit')
        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        molecule._off_molecule = openforcefield_toolkit.from_rdkit(molecule)

        molecule._build_rotamers()

        return molecule

    @staticmethod
    def from_openff(openff_molecule, rotamer_resolution=30,
                    exclude_terminal_rotamers=True, name='', tag='UNK',
                    connectivity_template=None, core_constraints=[],
                    allow_undefined_stereo=False):
        """
        It initializes and returns a peleffy Molecule representation
        from an OpenForceField molecular representation.

        Parameters
        ----------
        openff_molecule : an openforcefield.topology.Molecule object
            The OpenForceField's Molecule to use to initialize a peleffy
            Molecule object
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        name : str
            The molecule name
        tag : str
            The molecule tag. It must be a 3-character string
        connectivity_template : an rdkit.Chem.rdchem.Mol object
            A molecule represented with RDKit to use when assigning the
            connectivity of this Molecule object
        core_constraints : list[int or str]
            It defines the list of atoms to constrain in the core, thus,
            the core will be forced to contain them. Atoms can be specified
            through integers that match the atom index or strings that
            match with the atom PDB name
        allow_undefined_stereo : bool
            Whether to allow a molecule with undefined stereochemistry
            to be defined or try to assign the stereochemistry and
            raise a complaint if not possible. Default is False

        Returns
        -------
        molecule : an peleffy.topology.Molecule
            The resulting peleffy's Molecule object

        Examples
        --------

        Load a molecule from an RDKit molecular representation

        >>> from rdkit import Chem

        >>> rdkit_molecule = Chem.MolFromPDBFile(pdb_path)

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule.from_rdkit(rdkit_molecule)

        """
        if name == '':
            name = openff_molecule.name

        molecule = Molecule(
            rotamer_resolution=30,
            exclude_terminal_rotamers=exclude_terminal_rotamers,
            name=name, tag=tag,
            connectivity_template=connectivity_template,
            core_constraints=core_constraints,
            allow_undefined_stereo=allow_undefined_stereo)

        logger = Logger()

        logger.info(' - Initializing molecule from an OpenFF '
                    + 'molecular representation')
        molecule._initialize()
        molecule._off_molecule = openff_molecule

        logger.info('   - Generating RDKit molecular representation with '
                    + 'the Open Force Field Toolkit')
        molecule._rdkit_molecule = openff_molecule.to_rdkit()

        molecule._build_rotamers()

        return molecule

    @property
    def rotamer_resolution(self):
        """
        The resolution to be used when generating the rotamers library.

        Returns
        -------
        rotamer_resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space
        """
        return self._rotamer_resolution

    @property
    def exclude_terminal_rotamers(self):
        """
        The behavior when handling terminal rotamers when generating the
        rotamers library.

        Returns
        -------
        exclude_terminal_rotamers : bool
            Whether to exclude terminal rotamers when generating the
            rotamers library  or not
        """
        return self._exclude_terminal_rotamers

    @property
    def connectivity_template(self):
        """
        The template containing the correct connectivity for this Molecule
        object.

        Returns
        -------
        connectivity_template : an rdkit.Chem.rdchem.Mol object
            A molecule represented with RDKit to use when assigning the
            connectivity of this Molecule object
        """
        return self._connectivity_template

    @property
    def allow_undefined_stereo(self):
        """
        Whether to allow a molecule with undefined stereochemistry
        to be defined or try to assign the stereochemistry and
        raise a complaint if not possible

        Returns
        -------
        allow_undefined_stereo : bool
            The current configuration towards the stereochemistry
            behaviour of this molecule
        """
        return self._allow_undefined_stereo

    @property
    def fix_pdb(self):
        """
        It activates or deactivates the PDB fixer prior parsing the PDB
        file.

        Returns
        -------
        fix_pdb : bool
            The current activation state of the PDB fixer
        """
        return self._fix_pdb

    @property
    def core_constraints(self):
        """
        The list of indices or PDB names of the atoms to constraint to
        the core when building the rotamers.

        Returns
        -------
        core_constraints : list[int or str]
            The list of indices or PDB names of the atoms to constrain
        """
        return self._core_constraints

    @property
    def off_molecule(self):
        """
        The OpenForceField's molecule instance linked to the molecule.

        Returns
        -------
        off_molecule : an openforcefield.topology.Molecule
            The OpenForceField's molecule linked to this Molecule object
        """
        return self._off_molecule

    @property
    def rdkit_molecule(self):
        """
        The RDKit's molecule instance linked to the molecule.

        Returns
        -------
        rdkit_molecule : an rdkit.Chem.rdchem.Mol object
            The RDKit's molecule linked to this Molecule object
        """
        return self._rdkit_molecule

    @property
    def rotamers(self):
        """
        The list of rotamers of the molecule.

        Returns
        -------
        rotamers : list[list]
            The list of rotamers grouped by the branch they belong to
        """
        return self._rotamers

    @property
    def name(self):
        """
        Molecule's name.

        Returns
        -------
        name : str
            The name of this Molecule object
        """
        return self._name

    @property
    def tag(self):
        """
        Molecule's tag.

        Returns
        -------
        tag : str
            The tag of this Molecule object
        """
        return self._tag

    @property
    def graph(self):
        """
        The topological graph of the molecule.

        Returns
        -------
        graph : an peleffy.topology.rotamer.MolecularGraph object
            The topological graph of this Molecule object.
        """
        return self._graph

    def _ipython_display_(self):
        """
        It returns a RDKit molecule with an embeded 2D representation.

        Returns
        -------
        representation_2D : a IPython display object
            It is displayable RDKit molecule with an embeded 2D
            representation
        """
        from IPython.display import display

        # Get 2D molecular representation
        rdkit_toolkit = RDKitToolkitWrapper()
        representation = rdkit_toolkit.get_2D_representation(self)

        # Get its image
        image = rdkit_toolkit.draw_molecule(representation)

        return display(image)
