"""
This module handles all classes and functions related with molecular
representations.
"""

from pathlib import Path

from .topology import Bond, Angle, OFFProper, OFFImproper
from .rotamer import MolecularGraph
from offpele.utils.toolkits import (RDKitToolkitWrapper,
                                    OpenForceFieldToolkitWrapper)
from offpele.charge import (Am1bccCalculator, GasteigerCalculator)


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

    def __init__(self, path=None):
        """
        It initializes a Molecule object.

        Parameters
        ----------
        path : str
            The path to a PDB with the molecule structure

        Examples
        --------

        Load a molecule from a PDB file and parameterize it

        >>> from offpele.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')
        >>> molecule.parameterize('openff_unconstrained-1.1.1.offxml')

        Generate the rotamer library of a molecule

        >>> from offpele.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')
        >>> molecule.parameterize('openff_unconstrained-1.1.1.offxml')
        >>> molecule.build_rotamer_library(resolution=30)
        >>> molecule.rotamer_library.to_file('MOL.rot.assign')

        """
        if isinstance(path, str):
            from pathlib import Path
            extension = Path(path).suffix
            extension = extension.strip('.')
            if extension == 'pdb':
                self._initialize_from_pdb(path)
            else:
                raise ValueError(
                    '{} is not a valid extension'.format(extension))
        else:
            self._initialize()

    def _initialize(self):
        """It initializes an empty molecule."""
        self._name = ''
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
        self._rotamer_library = None
        self._graph = None

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

        # Set molecule name according to PDB's residue name
        name = rdkit_toolkit.get_residue_name(self)
        self.set_name(name)

        openforcefield_toolkit = OpenForceFieldToolkitWrapper()

        self._off_molecule = openforcefield_toolkit.from_rdkit(self)

    def set_name(self, name):
        """
        It sets the name of the molecule.

        Parameters
        ----------
        name : str
            The name to set to the molecule
        """
        if isinstance(name, str) and len(name) > 2:
            name = name[0:3].upper()
            self._name = name

            if self.off_molecule:
                self.off_molecule.name = name

    def parameterize(self, forcefield, charges_method=None):
        """
        It parameterizes the molecule with a certain forcefield.

        Parameters
        ----------
        forcefield : str or openforcefield.typing.engines.smirnoff.ForceField
                     object
            The forcefield from which the parameters will be obtained
        charges_method : str
            The name of the charges method to employ
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

        self._build_atoms()

        self._build_bonds()

        self._build_angles()

        self._build_propers()

        self._build_impropers()

    def build_rotamer_library(self, resolution=30, n_rot_bonds_to_ignore=1):
        """
        It builds the rotamer library of a parameterized molecule.

        .. todo ::

            * Consider moving this to the rotamer module.

        Parameters
        ----------
        resolution : float
            The resolution in degrees to discretize the rotamer's
            conformational space. Default is 30
        n_rot_bonds_to_ignore : int
            The number of terminal rotatable bonds to ignore when
            building the rotamer library. Default is 1
        """
        self._assert_parameterized()

        print(' - Generating rotamer library')

        self._graph = MolecularGraph(self)

        self.graph.set_core()

        self.graph.set_parents()

        self._rotamer_library = self.graph.build_rotamer_library(
            resolution, n_rot_bonds_to_ignore)

    def plot_rotamer_graph(self):
        """It plots the rotamer graph in screen."""
        self._assert_parameterized()

        try:
            from rdkit import Chem
        except ImportError:
            raise Exception('RDKit Python API not found')

        # Find rotatable bond ids as in Lipinski module in RDKit
        # https://github.com/rdkit/rdkit/blob/1bf6ef3d65f5c7b06b56862b3fb9116a3839b229/rdkit/Chem/Lipinski.py#L47
        rot_bonds_atom_ids = self._rdkit_molecule.GetSubstructMatches(
            Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]'))

        graph = self._compute_rotamer_graph(rot_bonds_atom_ids)

        rot_edges = [(u, v) for (u, v, d) in graph.edges(data=True)
                     if d['weight'] == 1]
        nrot_edges = [(u, v) for (u, v, d) in graph.edges(data=True)
                      if d['weight'] == 0]

        import networkx as nx

        pos = nx.circular_layout(graph)

        nx.draw_networkx_nodes(graph, pos, node_size=400)
        nx.draw_networkx_edges(graph, pos, edgelist=rot_edges,
                               width=4)
        nx.draw_networkx_edges(graph, pos, edgelist=nrot_edges,
                               width=4, alpha=0.5, edge_color='b',
                               style='dashed')
        nx.draw_networkx_labels(graph, pos, font_size=10,
                                font_family='sans-serif')

        import matplotlib.pyplot as plt
        plt.axis('off')
        plt.show()

    def _assert_parameterized(self):
        """
        It checks that the molecule has been parameterized, raises an
        AssertionError otherwise.
        """
        try:
            assert self.off_molecule is not None
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

    def _build_atoms(self):
        """It builds the atoms of the molecule."""
        # PELE needs underscores instead of whitespaces
        pdb_atom_names = {(i, ): name.replace(' ', '_',)
                          for i, name in enumerate(self.get_pdb_atom_names())}

        # TODO should we assign an OPLS type? How can we do this with OFF?
        OPLS_types = {i: None
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
            bond = Bond(index=index, atom1_idx=atom1_idx, atom2_idx=atom2_idx,
                        spring_constant=ks[atom_indexes],
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
            angle = Angle(index=index, atom1_idx=atom1_idx,
                          atom2_idx=atom2_idx, atom3_idx=atom3_idx,
                          spring_constant=ks[atom_indexes],
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

                if period and phase and k and idivf:
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

                if period and phase and k and idivf:
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

    def get_pdb_atom_names(self):
        """
        It returns the PDB atom names of all the atoms in the molecule.

        Returns
        -------
        pdb_atom_names : str
            The PDB atom names of all the atoms in this Molecule object
        """
        self._assert_parameterized()

        pdb_atom_names = list()

        for atom in self.rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            pdb_atom_names.append(pdb_info.GetName())

        return pdb_atom_names

    def to_impact(self, path):
        """
        .. todo ::

            * We still need to implement this
        """
        pass

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
    def rotamer_library(self):
        """
        The rotamer library of the molecule.

        Returns
        -------
        rotamer_library : an offpele.topology.rotamer.RotamerLibrary object
            The rotamer library of this Molecule object
        """
        return self._rotamer_library

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
