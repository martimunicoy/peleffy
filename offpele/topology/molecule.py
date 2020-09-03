"""
This module handles all classes and functions related with molecular
representations.
"""

from pathlib import Path

from .topology import Bond, Angle, OFFProper, OFFImproper
from .rotamer import MolecularGraph
from offpele.utils.toolkits import (RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper,
                                    SchrodingerToolkitWrapper)
from offpele.charge import (Am1bccCalculator, GasteigerCalculator,
                            OPLSChargeCalculator)


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
        parent : offpele.topology.Atom
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
        parent : offpele.topology.Atom
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
        parent : offpele.topology.Atom
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
                 exclude_terminal_rotamers=True, name='', tag='UNK'):
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

        Examples
        --------

        Load a molecule from a PDB file and parameterize it with Open
        Force Field

        >>> from offpele.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')
        >>> molecule.parameterize('openff_unconstrained-1.1.1.offxml')

        Load a molecule using a SMILES tag and parameterize it with Open
        Force Field

        >>> from offpele.topology import Molecule

        >>> molecule = Molecule(smiles='Cc1ccccc1')
        >>> molecule.parameterize('openff_unconstrained-1.2.0.offxml')

        """
        self._name = name
        self._tag = tag
        self._rotamer_resolution = rotamer_resolution
        self._exclude_terminal_rotamers = exclude_terminal_rotamers

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
        self._parameterized = False
        self._OPLS_included = False
        self._OPLS_parameters = None

    def _initialize_from_pdb(self, path):
        """
        It initializes a molecule with the molecule structure read from
        a PDB file.

        Parameters
        ----------
        path : str
            The path to a PDB with the molecule structure
        """
        self._initialize()
        print(' - Loading molecule from RDKit')

        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_pdb(path)

        # RDKit must generate stereochemistry specifically from 3D coords
        rdkit_toolkit.assign_stereochemistry_from_3D(self)

        # Set molecule name according to PDB name
        if self.name == '':
            from pathlib import Path
            name = Path(path).stem
            self.set_name(name)

        # Set molecule tag according to PDB's residue name
        if self.tag == 'UNK':
            tag = rdkit_toolkit.get_residue_name(self)
            self.set_tag(tag)

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
        self._initialize()
        print(' - Constructing molecule from a SMILES tag with RDKit')

        rdkit_toolkit = RDKitToolkitWrapper()
        self._rdkit_molecule = rdkit_toolkit.from_smiles(smiles)

        # TODO not sure if stereochemistry assignment from 3D is still necessary
        # RDKit must generate stereochemistry specifically from 3D coords
        # rdkit_toolkit.assign_stereochemistry_from_3D(self)

        # Set molecule name according to the SMILES tag
        if self.name == '':
            self.set_name(smiles)

        openforcefield_toolkit = OpenForceFieldToolkitWrapper()

        self._off_molecule = openforcefield_toolkit.from_rdkit(self)

    def _build_rotamers(self):
        """It builds the rotamers of the molecule."""
        if self.off_molecule and self.rdkit_molecule:
            print(' - Generating rotamer library')

            self._graph = MolecularGraph(self)
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

    def parameterize(self, forcefield, charges_method=None,
                     use_OPLS_nonbonding_params=False,
                     use_OPLS_bonds_and_angles=False):
        """
        It parameterizes the molecule with a certain forcefield.

        Parameters
        ----------
        forcefield : str or openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield from which the parameters will be obtained
        charges_method : str
            The name of the charges method to employ. One of
            ['gasteiger', 'am1bcc', 'OPLS']. If None, 'am1bcc' will be used
        use_OPLS_nonbonding_params : bool
            Whether to use Open Force Field or OPLS to obtain the
            nonbonding parameters. Please, note that this option is only
            available if a valid Schrodinger installation is found in the
            current machine. Default is False
        use_OPLS_bonds_and_angles : bool
            Whether to use OPLS to obtain the bond and angle parameters
            or not. Please, note that this option is only
            available if a valid Schrodinger installation is found in the
            current machine. Default is False
        """

        if not self.off_molecule or not self.rdkit_molecule:
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

        charges_calculator = self._get_charges_calculator(charges_method)

        print(' - Computing partial charges with '
              + '{}'.format(charges_calculator.name))
        self._assign_charges(charges_calculator)

        self._clean_lists()

        self._build_atoms()

        self._build_bonds()

        self._build_angles()

        self._build_propers()

        self._build_impropers()

        self._parameterized = True

        self.graph.set_core()

        self.graph.set_parents()

        self._OPLS_included = False

        if use_OPLS_nonbonding_params:
            self.add_OPLS_nonbonding_params()
        if use_OPLS_bonds_and_angles:
            self.add_OPLS_bonds_and_angles()

    def assert_parameterized(self):
        """
        It checks that the molecule has been parameterized, raises an
        AssertionError otherwise.
        """
        try:
            assert self.parameterized
        except AssertionError:
            raise Exception('Molecule not parameterized')

    def _get_charges_calculator(self, charges_method):
        """
        It returns the charges method that matches with the name that is
        supplied.

        Parameters
        ----------
        charges_method : str
            The name of the charges method to employ.
            One of ['gasteiger', 'am1bcc']. If None, 'am1bcc' will be used

        Returns
        -------
        charge_calculator : An offpele.topology.charges.PartialChargesCalculator
            object
            The charge calculation method that will be employed to calculate
            partial charges

        Raises
        ------
        Exception if the requested charge method is unknown
        """

        if charges_method == 'am1bcc' or charges_method is None:
            return Am1bccCalculator(self)

        elif charges_method == 'gasteiger':
            return GasteigerCalculator(self)

        elif charges_method == 'OPLS':
            return OPLSChargeCalculator(self)

        else:
            raise Exception('Charges method \'{}\' '.format(charges_method)
                            + 'is unknown')

    def _assign_charges(self, method):
        """It computes the partial charges using the charge calculation
        method that is supplied and assings them to the molecule..

        Parameters
        ----------
        method : An offpele.topology.charges.PartialChargesCalculator
            object
            The charge calculation method that will be employed to calculate
            partial charges
        """

        charges = method.get_partial_charges()

        self.off_molecule.partial_charges = charges

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
        # PELE needs underscores instead of whitespaces
        pdb_atom_names = {(i, ): name.replace(' ', '_',)
                          for i, name in enumerate(self.get_pdb_atom_names())}

        OPLS_types = {i: 'OFFT'
                      for i in self.parameters.get_vdW_parameters().keys()}

        # TODO Which is the purpose of unknown value? Is it important?
        unknowns = {i: None
                    for i in self.parameters.get_vdW_parameters().keys()}

        # TODO Create z-matrix from 3D coordinates

        coords = RDKitToolkitWrapper().get_coordinates(self)

        sigmas = self.parameters.get_vdW_sigmas()

        if all([sigma is None for sigma in sigmas.values()]):
            sigmas = self.parameters.get_vdW_sigmas_from_rmin_halves()

        epsilons = self.parameters.get_vdW_epsilons()

        # TODO Find a way to assign implicit solvent parameters to atoms with OFF
        born_radii = {i: None
                      for i in self.parameters.get_vdW_parameters().keys()}

        # TODO Doublecheck later this relation
        SASA_radii = {i: j / 2.0 for i, j in sigmas.items()}

        # TODO Find a way to assign implicit solvent parameters to atoms with OFF
        nonpolar_gammas = {i: None for i in
                           self.parameters.get_vdW_parameters().keys()}
        nonpolar_alphas = {i: None for i in
                           self.parameters.get_vdW_parameters().keys()}

        for index in self.parameters.get_vdW_parameters().keys():
            assert len(index) == 1, 'Index should be a tupple of length 1'
            atom = Atom(index=int(index[0]),
                        PDB_name=pdb_atom_names[index],
                        OPLS_type=OPLS_types[index],
                        unknown=unknowns[index],
                        x=coords[index][0],
                        y=coords[index][1],
                        z=coords[index][2],
                        sigma=sigmas[index],
                        epsilon=epsilons[index],
                        charge=self.off_molecule.partial_charges[index],
                        born_radius=born_radii[index],
                        SASA_radius=SASA_radii[index],
                        nonpolar_gamma=nonpolar_gammas[index],
                        nonpolar_alpha=nonpolar_alphas[index])
            self._add_atom(atom)

    def _add_atom(self, atom):
        """
        It adds an atom to the molecule's list of atoms.

        Parameters
        ----------
        atom : an offpele.topology.Atom
            The Atom to add
        """
        self._atoms.append(atom)

    def _build_bonds(self):
        """It builds the bonds of the molecule."""
        lengths = self.parameters.get_bond_lengths()

        ks = self.parameters.get_bond_ks()

        bond_indexes = self.parameters.get_bond_parameters().keys()

        for index, atom_indexes in enumerate(bond_indexes):
            (atom1_idx, atom2_idx) = atom_indexes
            # PELE works with half of the OFF's spring
            k = ks[atom_indexes] / 2.0
            bond = Bond(index=index, atom1_idx=atom1_idx, atom2_idx=atom2_idx,
                        spring_constant=k,
                        eq_dist=lengths[atom_indexes])
            self._add_bond(bond)

    def _add_bond(self, bond):
        """
        It adds a bond to the molecule's list of bonds.

        Parameters
        ----------
        bond : an offpele.topology.Bond
            The Bond to add
        """
        self._bonds.append(bond)

    def _build_angles(self):
        """It builds the angles of the molecule."""
        angles = self.parameters.get_angle_angles()

        ks = self.parameters.get_angle_ks()

        angle_indexes = self.parameters.get_angle_parameters().keys()

        for index, atom_indexes in enumerate(angle_indexes):
            atom1_idx, atom2_idx, atom3_idx = atom_indexes
            # PELE works with half of the OFF's spring
            k = ks[atom_indexes] / 2.0
            angle = Angle(index=index, atom1_idx=atom1_idx,
                          atom2_idx=atom2_idx, atom3_idx=atom3_idx,
                          spring_constant=k,
                          eq_angle=angles[atom_indexes])
            self._add_angle(angle)

    def _add_angle(self, angle):
        """
        It adds an angle to the molecule's list of angles.

        Parameters
        ----------
        angle : an offpele.topology.Angle
            The Angle to add
        """
        self._angles.append(angle)

    def _build_propers(self):
        """It builds the propers of the molecule."""
        periodicities = self.parameters.get_dihedral_periodicities()
        phases = self.parameters.get_dihedral_phases()
        ks = self.parameters.get_dihedral_ks()
        idivfs = self.parameters.get_dihedral_idivfs()

        # TODO in which situation these dicts are supposed to be None?
        if periodicities is None or phases is None or ks is None:
            return

        # idivf is a optional parameter in OpenForceField
        if len(idivfs) == 0:
            for period_by_index in periodicities:
                idivfs.append(dict(zip(
                    period_by_index.keys(),
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

                if (period is not None and phase is not None
                        and k is not None and idivf is not None):
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

        self._handle_excluded_propers()

    def _add_proper(self, proper):
        """
        It adds a proper dihedral to the molecule's list of propers.

        Parameters
        ----------
        proper : an offpele.topology.Proper
            The Proper to add
        """
        self._propers.append(proper)

    def _add_OFF_proper(self, proper):
        """
        It adds a proper dihedral to the molecule's list of OFF propers.

        Parameters
        ----------
        proper : an offpele.topology.OFFProper
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
        periodicities = self.parameters.get_improper_periodicities()
        phases = self.parameters.get_improper_phases()
        ks = self.parameters.get_improper_ks()
        idivfs = self.parameters.get_improper_idivfs()

        # TODO in which situation these dicts are supposed to be None?
        if periodicities is None or phases is None or ks is None:
            return

        # idivf is a optional parameter in OpenForceField
        if len(idivfs) == 0:
            for period_by_index in periodicities:
                idivfs.append(dict(zip(
                    period_by_index.keys(),
                    [1, ] * len(period_by_index.keys()))))

        assert len(periodicities) == len(phases) and \
            len(periodicities) == len(ks) and \
            len(periodicities) == len(idivfs), 'Unconsistent set of ' \
            'OpenForceField\'s improper parameters. They all should ' \
            'have equal lengths'

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

                if (period is not None and phase is not None
                        and k is not None and idivf is not None):
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
        """
        It adds an improper dihedral to the molecule's list of impropers.

        Parameters
        ----------
        improper : an offpele.topology.Improper
            The Improper to add
        """
        self._impropers.append(improper)

    def _add_OFF_improper(self, improper):
        """
        It adds an improper dihedral to the molecule's list of OFF impropers.

        Parameters
        ----------
        improper : an offpele.topology.OFFImproper
            The OFFImproper to add
        """
        self._OFF_impropers.append(improper)

    def add_OPLS_nonbonding_params(self):
        """
        It adds OPLS' nonbonding parameters to the molecule. Please, note
        that OPLS' partial charges are not set in this function.
        Instead, they are assigned in the Molecule's parameterize()
        function when 'OPLS' is chosen as the 'charges_method'.
        """

        self.assert_parameterized()

        OPLS_params = self.get_OPLS_parameters()

        for atom, atom_type, sigma, epsilon, SGB_radius, \
            vdW_radius, gamma, alpha in zip(self.atoms,
                                            OPLS_params['atom_types'],
                                            OPLS_params['sigmas'],
                                            OPLS_params['epsilons'],
                                            OPLS_params['SGB_radii'],
                                            OPLS_params['vdW_radii'],
                                            OPLS_params['gammas'],
                                            OPLS_params['alphas']):
            atom.set_OPLS_type(atom_type)
            atom.set_sigma(sigma)
            atom.set_epsilon(epsilon)
            atom.set_born_radius(SGB_radius)
            atom.set_SASA_radius(vdW_radius)
            atom.set_nonpolar_gamma(gamma)
            atom.set_nonpolar_alpha(alpha)

        self._OPLS_included = True

    def add_OPLS_bonds_and_angles(self):
        """
        It adds OPLS' bond and angle parameters to the molecule.
        """

        self.assert_parameterized()

        OPLS_params = self.get_OPLS_parameters()

        self._bonds = list()
        for index, bond in enumerate(OPLS_params['bonds']):
            atom1 = self.atoms[bond['atom1_idx']]
            atom2 = self.atoms[bond['atom2_idx']]
            OPLS_bond = Bond(index=index,
                             atom1_idx=atom1.index,
                             atom2_idx=atom2.index,
                             spring_constant=bond['spring_constant'],
                             eq_dist=bond['eq_dist'])
            self._add_bond(OPLS_bond)

        self._angles = list()
        for index, angle in enumerate(OPLS_params['angles']):
            atom1 = self.atoms[angle['atom1_idx']]
            atom2 = self.atoms[angle['atom2_idx']]
            atom3 = self.atoms[angle['atom3_idx']]
            angle = Angle(index=index,
                          atom1_idx=atom1.index,
                          atom2_idx=atom2.index,
                          atom3_idx=atom3.index,
                          spring_constant=angle['spring_constant'],
                          eq_angle=angle['eq_angle'])
            self._add_angle(angle)

        self._OPLS_included = True

    def get_pdb_atom_names(self):
        """
        It returns the PDB atom names of all the atoms in the molecule.

        Returns
        -------
        pdb_atom_names : str
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

    def get_OPLS_parameters(self):
        """
        It returns the OPLS parameters of the molecule. It first looks
        if they have already been calculated and returns them if found.
        Otherwise, it uses the SchrodingerToolkitWrapper to generate
        them.

        Returns
        -------
        OPLS_parameters : a SchrodingerToolkitWrapper.OPLSParameters object
            The set of lists of parameters grouped by parameter type.
            Thus, the dictionary has the following keys: atom_names,
            atom_types, charges, sigmas, epsilons, SGB_radii, vdW_radii,
            gammas, and alphas
        """

        if self._OPLS_parameters is None:
            schrodinger_toolkit = SchrodingerToolkitWrapper()

            self._OPLS_parameters = \
                schrodinger_toolkit.get_OPLS_parameters(self)

        return self._OPLS_parameters

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
        The forcefield employed to parameterize the molecule.

        Returns
        -------
        forcefield : an openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield employed to parameterize this Molecule object
        """
        return self._forcefield

    @property
    def atoms(self):
        """
        The list of atoms of the molecule.

        Returns
        -------
        atoms : list[offpele.topology.molecule.Atom]
            The list of atoms of this Molecule object.
        """
        return self._atoms

    @property
    def bonds(self):
        """
        The list of bonds of the molecule.

        Returns
        -------
        bonds : list[offpele.topology.Bond]
            The list of bonds of this Molecule object.
        """
        return self._bonds

    @property
    def angles(self):
        """
        The list of angles of the molecule.

        Returns
        -------
        angles : list[offpele.topology.Angle]
            The list of angles of this Molecule object.
        """
        return self._angles

    @property
    def propers(self):
        """
        The list of propers of the molecule.

        Returns
        -------
        propers : list[offpele.topology.Proper]
            The list of propers of this Molecule object.
        """
        return self._propers

    @property
    def impropers(self):
        """
        The list of impropers of the molecule.

        Returns
        -------
        impropers : list[offpele.topology.Improper]
            The list of impropers of this Molecule object.
        """
        return self._impropers

    @property
    def graph(self):
        """
        The topological graph of the molecule.

        Returns
        -------
        graph : an offpele.topology.rotamer.MolecularGraph object
            The topological graph of this Molecule object.
        """
        return self._graph

    @property
    def parameterized(self):
        """
        Whether the molecule has been parameterized with the Open Force Field
        toolkit or not.

        Returns
        -------
        parameterized : bool
            The parameterization status
        """
        return self._parameterized

    @property
    def OPLS_included(self):
        """
        Whether the molecule has been parameterized with OPLS in combination
        with the Open Force Field toolkit or not.

        Returns
        -------
        OPLS_included : bool
            The OPLS combination status
        """
        return self._OPLS_included

    @property
    def OPLS_parameters(self):
        """
        The OPLS parameters of the molecule.

        Returns
        -------
        OPLS_parameters : a SchrodingerToolkitWrapper.OPLSParameters object
            The set of lists of parameters grouped by parameter type.
            Thus, the dictionary has the following keys: atom_names,
            atom_types, charges, sigmas, epsilons, SGB_radii, vdW_radii,
            gammas, and alphas
        """
        return self.get_OPLS_parameters()

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

        rdkit_toolkit = RDKitToolkitWrapper()
        return display(rdkit_toolkit.get_2D_representation(self))
