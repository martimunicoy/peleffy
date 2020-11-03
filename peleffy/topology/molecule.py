"""
This module handles all classes and functions related with molecular
representations.
"""


from .topology import Bond, Angle, OFFProper, OFFImproper
from .rotamer import MolecularGraph, MolecularGraphWithConstrainedCore
from peleffy.utils.toolkits import (RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper)
from peleffy.forcefield import ForceFieldSelector
from peleffy.charge import (Am1bccCalculator, GasteigerCalculator,
                            OPLSChargeCalculator)
from peleffy.utils import Logger


class Atom(object):
    """
    It represents the properties of an atom.
    """

    def __init__(self, index=-1, core=None, OPLS_type=None, PDB_name=None,
                 unknown=None, x=None, y=None, z=None, sigma=None,
                 epsilon=None, charge=None, born_radius=None, SASA_radius=None,
                 nonpolar_gamma=None, nonpolar_alpha=None, parent=None):
        """
        It initializes an Atom object.

        Parameters
        ----------
        index : int
            The index of the atom
        core : bool
            Whether atom is in the core or in a branch
        OPLS_type : str
            The OPLS type of the atom
        PDB_name : str
            The PDB name of the atom
        unknown : int
            The unknown value of the atom
        x : float
            The x coordinate of the atom
        y : float
            The y coordinate of the atom
        z : float
            The z coordinate of the atom
        sigma : simtk.unit.Quantity
            The sigma parameter of the atom
        epsilon : simtk.unit.Quantity
            The epsilon parameter of the atom
        charge : simtk.unit.Quantity
            The partial charge parameter of the atom
        born_radius : simtk.unit.Quantity
            The SGB born parameter radius of the atom
        SASA_radius : simtk.unit.Quantity
            The SASA radius parameter of the atom
        nonpolar_gamma : simtk.unit.Quantity
            The nonpolar gamma parameter of the atom
        nonpolar_alpha : simtk.unit.Quantity
            The nonpolar alpha parameter of the atom
        parent : peleffy.topology.Atom
            The parent of the atom
        """
        self._index = index
        self._core = core
        self._OPLS_type = OPLS_type
        self._PDB_name = PDB_name
        self._unknown = unknown
        self._x = x
        self._y = y
        self._z = z
        self._sigma = sigma
        self._epsilon = epsilon
        self._charge = charge
        self._born_radius = born_radius  # Rad. Non Polar SGB
        self._SASA_radius = SASA_radius  # Rad. Non Polar Type
        self._nonpolar_gamma = nonpolar_gamma  # SGB Non Polar gamma
        self._nonpolar_alpha = nonpolar_alpha  # SGB Non Polar type
        self._parent = parent

    def set_index(self, index):
        """
        It sets the index of the atom.

        Parameters
        ----------
        index : int
            The index of this Atom object
        """
        self._index = index

    def set_as_core(self):
        """It sets the atom as core"""
        self._core = True

    def set_as_branch(self):
        """It sets the atom as branch"""
        self._core = False

    def set_parent(self, parent):
        """
        It sets the parent of the atom.

        Parameters
        ----------
        parent : peleffy.topology.Atom
            The parent of the atom
        """
        self._parent = parent

    def set_coords(self, coords):
        """
        It sets the coordinates of the atom.

        Parameters
        ----------
        coords : list
            The coordinates array to set to this Atom object
        """
        assert len(coords) == 3, '3D array is expected'

        self._x, self._y, self._z = coords

    def set_OPLS_type(self, OPLS_type):
        """
        It sets the OPLS type of the atom.

        Parameters
        ----------
        OPLS_type : str
            The OPLS type to set to this Atom object
        """
        self._OPLS_type = OPLS_type

    def set_sigma(self, sigma):
        """
        It sets the sigma of the atom.

        Parameters
        ----------
        sigma : float
            The sigma to set to this Atom object
        """
        self._sigma = sigma

    def set_epsilon(self, epsilon):
        """
        It sets the epsilon of the atom.

        Parameters
        ----------
        epsilon : float
            The epsilon to set to this Atom object
        """
        self._epsilon = epsilon

    def set_charge(self, charge):
        """
        It sets the charge of the atom.

        Parameters
        ----------
        charge : float
            The charge to set to this Atom object
        """
        self._charge = charge

    def set_born_radius(self, born_radius):
        """
        It sets the Born radius of the atom.

        Parameters
        ----------
        born_radius : float
            The Born radius to set to this Atom object
        """
        self._born_radius = born_radius

    def set_SASA_radius(self, SASA_radius):
        """
        It sets the SASA radius of the atom.

        Parameters
        ----------
        SASA_radius : float
            The SASA radius to set to this Atom object
        """
        self._SASA_radius = SASA_radius

    def set_nonpolar_gamma(self, nonpolar_gamma):
        """
        It sets the nonpolar gamma of the atom.

        Parameters
        ----------
        nonpolar_gamma : float
            The nonpolar gamma to set to this Atom object
        """
        self._nonpolar_gamma = nonpolar_gamma

    def set_nonpolar_alpha(self, nonpolar_alpha):
        """
        It sets the nonpolar alpha of the atom.

        Parameters
        ----------
        nonpolar_alpha : float
            The nonpolar alpha to set to this Atom object
        """
        self._nonpolar_alpha = nonpolar_alpha

    @property
    def index(self):
        """
        Atom's index.

        Returns
        -------
        index : int
            The index of this Atom object
        """
        return self._index

    @property
    def core(self):
        """
        Atom's core position.

        Returns
        -------
        core : bool
            Whether this Atom object is in the core or in a branch
        """
        return self._core

    @property
    def OPLS_type(self):
        """
        Atom's OPLS type.

        .. todo ::

           * Consider removing any reference to OPLS, if possible
             Otherwise, use SMIRks to find the best match

        Returns
        -------
        OPLS_type : str
            The OLPS type of this Atom object
        """
        return self._OPLS_type

    @property
    def PDB_name(self):
        """
        Atom's PDB name.

        .. todo ::

           * Consider removing any reference to OPLS, if possible
             Otherwise, use SMIRks to find the best match

        Returns
        -------
        PDB_name : str
            The PDB name of this Atom object
        """
        return self._PDB_name

    @property
    def unknown(self):
        """
        Atom's unknown int.

        .. todo ::

           * Review the actual purpose of this attribute in PELE

        Returns
        -------
        unknown : int
            The unknown int of this Atom object
        """
        return self._unknown

    @property
    def x(self):
        """
        Atom's x coordinate.

        Returns
        -------
        x : float
            The x coordinate of this Atom object
        """
        return self._x

    @property
    def y(self):
        """
        Atom's y coordinate.

        Returns
        -------
        y : float
            The y coordinate of this Atom object
        """
        return self._y

    @property
    def z(self):
        """
        Atom's z coordinate.

        Returns
        -------
        z : float
            The z coordinate of this Atom object
        """
        return self._z

    @property
    def sigma(self):
        """
        Atom's sigma.

        Returns
        -------
        sigma : simtk.unit.Quantity
            The sigma parameter of this Atom object
        """
        return self._sigma

    @property
    def epsilon(self):
        """
        Atom's epsilon.

        Returns
        -------
        epsilon : simtk.unit.Quantity
            The epsilon parameter of this Atom object
        """
        return self._epsilon

    @property
    def charge(self):
        """
        Atom's charge.

        Returns
        -------
        charge : simtk.unit.Quantity
            The charge parameter of this Atom object
        """
        return self._charge

    @property
    def born_radius(self):
        """
        Atom's born radius.

        Returns
        -------
        born_radius : simtk.unit.Quantity
            The SGB Born radius parameter of this Atom object
        """
        return self._born_radius

    @property
    def SASA_radius(self):
        """
        Atom's SASA radius.

        Returns
        -------
        SASA_radius : simtk.unit.Quantity
            The SASA radius parameter of this Atom object
        """
        return self._SASA_radius

    @property
    def nonpolar_gamma(self):
        """
        Atom's nonpolar gamma.

        Returns
        -------
        nonpolar_gamma : simtk.unit.Quantity
            The nonpolar gamma parameter of this Atom object
        """
        return self._nonpolar_gamma

    @property
    def nonpolar_alpha(self):
        """
        Atom's nonpolar alpha.

        Returns
        -------
        nonpolar_alpha : simtk.unit.Quantity
            The nonpolar alpha parameter of this Atom object
        """
        return self._nonpolar_alpha

    @property
    def parent(self):
        """
        Atom's parent.

        Returns
        -------
        parent : simtk.unit.Quantity
            The nonpolar gamma parameter of this Atom object
        """
        return self._parent


class DummyAtom(Atom):
    """
    It represents a dummy atom.
    """

    def __init__(self, index=-1, PDB_name='DUMM', parent=None):
        """
        It initializes a DummyAtom object.

        Parameters
        ----------
        index : int
            The index of the atom
        PDB_name : str
            The PDB name of the atom
        parent : peleffy.topology.Atom
            The parent of the atom
        """
        if parent is None:
            parent = self
        super().__init__(index, False, None, PDB_name, None, None, None, None,
                         None, None, None, None, None, None, None, parent)


class Molecule(object):
    """
    It represent wraps up all the tools to parameterize a molecule with
    the OpenForceField toolkit for PELE.
    """

    def __init__(self, path=None, smiles=None, rotamer_resolution=30,
                 exclude_terminal_rotamers=True, name='', tag='UNK',
                 connectivity_template=None, core_constraints=[]):
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
        >>> molecule.parameterize('openff_unconstrained-1.2.0.offxml')

        Load a molecule usign a PDB file (without connectivity) and assign
        the missing connectivity from an RDKit template (e.g. obtained
        from qcportal and the Open Force Field Toolkit)

        >>> import qcportal as ptl
        >>> from openforcefield.topology import Molecule as OFFMolecule

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

        """
        self._name = name
        self._tag = tag
        self._rotamer_resolution = rotamer_resolution
        self._exclude_terminal_rotamers = exclude_terminal_rotamers
        self._connectivity_template = connectivity_template
        self._core_constraints = core_constraints

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

        else:
            self._initialize()

        self._build_rotamers()

    def _initialize(self):
        """It initializes an empty molecule."""
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
        self._rotamers = None
        self._graph = None

        from peleffy.forcefield.parameters import BaseParameterWrapper
        self._parameters = BaseParameterWrapper()

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

        logger.info('   - Loading molecule from RDKit')
        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_pdb(path)

        # Use RDKit template, if any, to assign the connectivity to
        # the current Molecule object
        if self.connectivity_template is not None:
            logger.info('   - Assigning connectivity from template')
            rdkit_toolkit.assign_connectivity_from_template(self)

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

    def set_forcefield(self, forcefield):
        """
        It sets the force field of the molecule.

        Parameters
        ----------
        forcefield : any _BaseForceField object
            The force field to employ to parameterize the molecule
        """
        from peleffy.forcefield import OpenForceField as ff1
        from peleffy.forcefield import OPLS2005ForceField as ff2
        from peleffy.forcefield import OpenFFOPLS2005ForceField as ff3
        from peleffy.forcefield.forcefield \
            import OpenForceField as ff4
        from peleffy.forcefield.forcefield \
            import OPLS2005ForceField as ff5
        from peleffy.forcefield.forcefield \
            import OpenFFOPLS2005ForceField as ff6

        if not isinstance(forcefield, (ff1, ff2, ff3, ff4, ff5, ff6)):
            raise TypeError('An invalid force field was supplied')

        self._forcefield = forcefield

    def parameterize(self, forcefield_name=None, charge_method=None,
                     force_parameterization=False):
        """
        It parameterizes the molecule with a certain forcefield.

        Parameters
        ----------
        forcefield_name : str
            The name of the forcefield from which the parameters will be
            obtained
        charge_method : str
            The name of the charge method to employ. One of
            ['gasteiger', 'am1bcc', 'OPLS']. If None, 'am1bcc' will be used
        force_parameterization : bool
            Whether to force a new parameterization instead of attempting
            to reuse parameters obtained in a previous parameterization,
            or not
        """

        if not self.off_molecule or not self.rdkit_molecule:
            raise Exception('OpenForceField molecule was not initialized '
                            + 'correctly')

        logger = Logger()
        logger.info(' - Loading forcefield')
        ff_selector = ForceFieldSelector()

        # Set forcefield and the corresponding parameters
        if forcefield_name is not None:
            forcefield = ff_selector.get_by_name(forcefield_name)
            self.set_forcefield(forcefield)
        else:
            if self.forcefield is None:
                raise ValueError('No force field has been set')
        self._parameters = self.forcefield.parameterize(self,
                                                        force_parameterization)

        # Initialize the charges calculator
        charge_calculator = self._get_charge_calculator(charge_method)

        logger.info(' - Computing partial charges with '
                    + '{}'.format(charge_calculator.name))
        self._assign_charges(charge_calculator)

        self._clean_lists()

        self._build_atoms()

        self._build_bonds()

        self._build_angles()

        self._build_propers()

        self._build_impropers()

        self.graph.set_core()

        self.graph.set_parents()

    def assert_parameterized(self):
        """
        It checks that the molecule has been parameterized, raises an
        AssertionError otherwise.
        """
        try:
            assert self.parameterized
        except AssertionError:
            raise Exception('Molecule not parameterized')

    def _get_charge_calculator(self, charge_method):
        """
        It returns the charge method that matches with the name that is
        supplied.

        .. todo ::

            * Move this function to an external charge calculator selector

        Parameters
        ----------
        charge_method : str
            The name of the charges method to employ.
            One of ['gasteiger', 'am1bcc']. If None, 'am1bcc' will be used

        Returns
        -------
        charge_calculator : An peleffy.topology.charges.PartialChargesCalculator
            object
            The charge calculation method that will be employed to calculate
            partial charges

        Raises
        ------
        ValueError
            If the requested charge method is unknown
        """

        if charge_method == 'am1bcc':
            return Am1bccCalculator(self)

        elif charge_method == 'gasteiger':
            return GasteigerCalculator(self)

        elif charge_method == 'OPLS':
            return OPLSChargeCalculator(self)

        elif charge_method is None:
            if self.forcefield.type == 'OPLS2005':
                return OPLSChargeCalculator(self)
            else:
                return Am1bccCalculator(self)

        else:
            raise ValueError('Charge method \'{}\' '.format(charge_method)
                             + 'is unknown')

    def _assign_charges(self, method):
        """It computes the partial charges using the charge calculation
        method that is supplied and assings them to the molecule..

        Parameters
        ----------
        method : An peleffy.topology.charges.PartialChargesCalculator
            object
            The charge calculation method that will be employed to calculate
            partial charges
        """

        charges = method.get_partial_charges()

        self.parameters['charges'] = charges

    def _clean_lists(self):
        """It cleans all the lists before parameterizing."""
        self._atoms = list()
        self._bonds = list()
        self._angles = list()
        self._propers = list()
        self._OFF_propers = list()
        self._impropers = list()
        self._OFF_impropers = list()

    def _build_atoms(self):
        """It builds the atoms of the molecule."""
        coords = RDKitToolkitWrapper().get_coordinates(self)

        for index, (atom_name, atom_type, sigma, epsilon, charge,
                    SGB_radius, vdW_radius, gamma, alpha) \
                in enumerate(self.parameters.atom_iterator):
            atom = Atom(index=index,
                        PDB_name=atom_name,
                        OPLS_type=atom_type,
                        x=coords[index][0],
                        y=coords[index][1],
                        z=coords[index][2],
                        sigma=sigma,
                        epsilon=epsilon,
                        charge=charge,
                        born_radius=SGB_radius,
                        SASA_radius=vdW_radius,
                        nonpolar_gamma=gamma,
                        nonpolar_alpha=alpha)
            self._add_atom(atom)

    def _add_atom(self, atom):
        """
        It adds an atom to the molecule's list of atoms.

        Parameters
        ----------
        atom : an peleffy.topology.Atom
            The Atom to add
        """
        self._atoms.append(atom)

    def _build_bonds(self):
        """It builds the bonds of the molecule."""
        for index, bond in enumerate(self.parameters['bonds']):
            bond = Bond(index=index,
                        atom1_idx=bond['atom1_idx'],
                        atom2_idx=bond['atom2_idx'],
                        spring_constant=bond['spring_constant'],
                        eq_dist=bond['eq_dist'])
            self._add_bond(bond)

    def _add_bond(self, bond):
        """
        It adds a bond to the molecule's list of bonds.

        Parameters
        ----------
        bond : an peleffy.topology.Bond
            The Bond to add
        """
        self._bonds.append(bond)

    def _build_angles(self):
        """It builds the angles of the molecule."""
        for index, angle in enumerate(self.parameters['angles']):
            angle = Angle(index=index,
                          atom1_idx=angle['atom1_idx'],
                          atom2_idx=angle['atom2_idx'],
                          atom3_idx=angle['atom3_idx'],
                          spring_constant=angle['spring_constant'],
                          eq_angle=angle['eq_angle'])
            self._add_angle(angle)

    def _add_angle(self, angle):
        """
        It adds an angle to the molecule's list of angles.

        Parameters
        ----------
        angle : an peleffy.topology.Angle
            The Angle to add
        """
        self._angles.append(angle)

    def _build_propers(self):
        """It builds the propers of the molecule."""
        for index, proper in enumerate(self.parameters['propers']):
            off_proper = OFFProper(atom1_idx=proper['atom1_idx'],
                                   atom2_idx=proper['atom2_idx'],
                                   atom3_idx=proper['atom3_idx'],
                                   atom4_idx=proper['atom4_idx'],
                                   periodicity=proper['periodicity'],
                                   phase=proper['phase'],
                                   k=proper['k'],
                                   idivf=proper['idivf'])

            PELE_proper = off_proper.to_PELE()
            self._add_proper(PELE_proper)
            self._add_OFF_proper(off_proper)

        self._handle_excluded_propers()

    def _add_proper(self, proper):
        """
        It adds a proper dihedral to the molecule's list of propers.

        Parameters
        ----------
        proper : an peleffy.topology.Proper
            The Proper to add
        """
        self._propers.append(proper)

    def _add_OFF_proper(self, proper):
        """
        It adds a proper dihedral to the molecule's list of OFF propers.

        Parameters
        ----------
        proper : an peleffy.topology.OFFProper
            The OFFProper to add
        """
        self._OFF_propers.append(proper)

    def _handle_excluded_propers(self):
        """
        It looks for those propers that define duplicated 1-4 relations
        and sets them to be ignored in PELE's 1-4 list.
        """
        for i, proper in enumerate(self.propers):
            atom1_idx = proper.atom1_idx
            atom4_idx = proper.atom4_idx
            for proper_to_compare in self.propers[0:i]:
                if proper == proper_to_compare:
                    continue

                if proper_to_compare.atom3_idx < 0:
                    continue

                # PELE already ignores 1-4 pair when the proper is exactly
                # the same
                if (proper.atom1_idx == proper_to_compare.atom1_idx
                        and proper.atom2_idx == proper_to_compare.atom2_idx
                        and proper.atom3_idx == proper_to_compare.atom3_idx
                        and proper.atom4_idx == proper_to_compare.atom4_idx):
                    continue

                atom1_idx_to_compare = proper_to_compare.atom1_idx
                atom4_idx_to_compare = proper_to_compare.atom4_idx
                if (atom1_idx == atom1_idx_to_compare
                        and atom4_idx == atom4_idx_to_compare):
                    proper.exclude_from_14_list()
                elif (atom1_idx == atom4_idx_to_compare
                        and atom4_idx == atom1_idx_to_compare):
                    proper.exclude_from_14_list()

            for angle_to_compare in self.angles:
                atom1_idx_to_compare = angle_to_compare.atom1_idx
                atom4_idx_to_compare = angle_to_compare.atom3_idx
                if (atom1_idx == atom1_idx_to_compare
                        and atom4_idx == atom4_idx_to_compare):
                    proper.exclude_from_14_list()
                elif (atom1_idx == atom4_idx_to_compare
                        and atom4_idx == atom1_idx_to_compare):
                    proper.exclude_from_14_list()

            for bond_to_compare in self.bonds:
                atom1_idx_to_compare = bond_to_compare.atom1_idx
                atom4_idx_to_compare = bond_to_compare.atom2_idx
                if (atom1_idx == atom1_idx_to_compare
                        and atom4_idx == atom4_idx_to_compare):
                    proper.exclude_from_14_list()
                elif (atom1_idx == atom4_idx_to_compare
                        and atom4_idx == atom1_idx_to_compare):
                    proper.exclude_from_14_list()

    def _build_impropers(self):
        """It builds the impropers of the molecule."""
        for index, improper in enumerate(self.parameters['impropers']):
            off_improper = OFFImproper(atom1_idx=improper['atom1_idx'],
                                       atom2_idx=improper['atom2_idx'],
                                       atom3_idx=improper['atom3_idx'],
                                       atom4_idx=improper['atom4_idx'],
                                       periodicity=improper['periodicity'],
                                       phase=improper['phase'],
                                       k=improper['k'],
                                       idivf=improper['idivf'])

            PELE_improper = off_improper.to_PELE()
            self._add_improper(PELE_improper)
            self._add_OFF_improper(off_improper)

    def _add_improper(self, improper):
        """
        It adds an improper dihedral to the molecule's list of impropers.

        Parameters
        ----------
        improper : an peleffy.topology.Improper
            The Improper to add
        """
        self._impropers.append(improper)

    def _add_OFF_improper(self, improper):
        """
        It adds an improper dihedral to the molecule's list of OFF impropers.

        Parameters
        ----------
        improper : an peleffy.topology.OFFImproper
            The OFFImproper to add
        """
        self._OFF_impropers.append(improper)

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

    def to_impact_file(self, path):
        """
        .. todo ::

            * We still need to implement this
        """
        # assert parameterized, then write impact
        pass

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
    def forcefield(self):
        """
        The name of the forcefield employed to parameterize the molecule.

        Returns
        -------
        forcefield : str
            The forcefield name
        """
        return self._forcefield

    @property
    def atoms(self):
        """
        The list of atoms of the molecule.

        Returns
        -------
        atoms : list[peleffy.topology.molecule.Atom]
            The list of atoms of this Molecule object.
        """
        return self._atoms

    @property
    def bonds(self):
        """
        The list of bonds of the molecule.

        Returns
        -------
        bonds : list[peleffy.topology.Bond]
            The list of bonds of this Molecule object.
        """
        return self._bonds

    @property
    def angles(self):
        """
        The list of angles of the molecule.

        Returns
        -------
        angles : list[peleffy.topology.Angle]
            The list of angles of this Molecule object.
        """
        return self._angles

    @property
    def propers(self):
        """
        The list of propers of the molecule.

        Returns
        -------
        propers : list[peleffy.topology.Proper]
            The list of propers of this Molecule object.
        """
        return self._propers

    @property
    def impropers(self):
        """
        The list of impropers of the molecule.

        Returns
        -------
        impropers : list[peleffy.topology.Improper]
            The list of impropers of this Molecule object.
        """
        return self._impropers

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

    @property
    def parameters(self):
        """
        It contains the parameter wrapper of the molecule. If the
        molecule has not been parameterized yet, it is set to None.

        Returns
        -------
        parameters : an BaseParameterWrapper object
            The parameter wrapper
        """
        return self._parameters

    @property
    def parameterized(self):
        """
        Whether the molecule has been parameterized or not.

        Returns
        -------
        parameterized : bool
            The parameterization status
        """
        return not self.parameters.is_empty()

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

        return display(representation)
