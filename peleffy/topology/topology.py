"""
This module handles the topological elements of force fields.
"""


from peleffy.topology.elements import (Atom, Bond, Angle,
                                       OFFProper, OFFImproper)
from peleffy.utils.toolkits import RDKitToolkitWrapper
from peleffy.utils import Logger


class Topology(object):
    """
    It represents the topology of a molecule.
    """

    def __init__(self, molecule, parameters):
        """
        It initializes a molecular topology representation, given a
        molecule and its parameters.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object whose topology will be built
        parameters : a BaseParameterWrapper
            The parameter wrapper belonging to the molecule

        Examples
        --------

        Load a molecule from a PDB file, parameterize it with a force
        field and generate its topology

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')
        >>> molecule.parameterize('openff_unconstrained-1.2.0.offxml')

        >>> from peleffy.forcefield import OpenForceField

        >>> openff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        >>> parameters = openff.parameterize(molecule)

        >>> from peleffy.topology import Topology

        >>> topology = Topology(molecule, parameters)

        Once a topology has been built, topological elements such as
        atoms, bonds, angles and dihedrals are easily accessible
        through the attributes below

        >>> topology.atoms

        >>> topology.bonds

        >>> topology.angles

        >>> topology.propers

        >>> topology.impropers

        """

        # Set topology attributes
        self._molecule = molecule
        self._parameters = parameters

        # Initialize and build the topology
        self._initialize()
        self._build()

    def _initialize(self):
        """It initializes all the lists before building the topology."""
        self._atoms = list()
        self._bonds = list()
        self._angles = list()
        self._propers = list()
        self._OFF_propers = list()
        self._impropers = list()
        self._OFF_impropers = list()

    def _build(self):
        """The topology builder."""

        # In case the molecule has not been initialized
        if (self.molecule.rdkit_molecule is None
                or len(list(self.parameters.atom_iterator)) == 0):
            logger = Logger()
            logger.warning('Warning: the input molecule has not been '
                           + ' initialized and its topology will be empty')
            return

        self._build_atoms()
        self._build_bonds()
        self._build_angles()
        self._build_propers()
        self._build_impropers()

    def _build_atoms(self):
        """It builds the atoms of the molecule."""

        from peleffy.utils import Logger
        logger = Logger()

        coords = RDKitToolkitWrapper().get_coordinates(self.molecule)

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
            self.add_atom(atom)

        for atom in self.atoms:
            if atom.index in self.molecule.graph.core_nodes:
                atom.set_as_core()
            else:
                atom.set_as_branch()

        # Start from an atom from the core
        absolute_parent = None
        for atom in self.atoms:
            if atom.core:
                absolute_parent = atom.index
                break
        else:
            logger.error('Error: no core atom found in molecule '
                         + '{}'.format(self.molecule.name))

        # Get parent indexes from the molecular graph
        parent_idxs = self.molecule.graph.get_parents(absolute_parent)

        # Assert parent_idxs has right length
        if len(parent_idxs) != len(self.atoms):
            logger.error('Error: no core atom found in molecule '
                         + '{}'.format(self.molecule.name))

        for atom in self.atoms:
            parent_idx = parent_idxs[atom.index]
            if parent_idx is not None:
                atom.set_parent(self.atoms[parent_idx])

    def _build_bonds(self):
        """It builds the bonds of the molecule."""
        for index, bond in enumerate(self.parameters['bonds']):
            bond = Bond(index=index,
                        atom1_idx=bond['atom1_idx'],
                        atom2_idx=bond['atom2_idx'],
                        spring_constant=bond['spring_constant'],
                        eq_dist=bond['eq_dist'])
            self.add_bond(bond)

    def _build_angles(self):
        """It builds the angles of the molecule."""
        for index, angle in enumerate(self.parameters['angles']):
            angle = Angle(index=index,
                          atom1_idx=angle['atom1_idx'],
                          atom2_idx=angle['atom2_idx'],
                          atom3_idx=angle['atom3_idx'],
                          spring_constant=angle['spring_constant'],
                          eq_angle=angle['eq_angle'])
            self.add_angle(angle)

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
            self.add_proper(PELE_proper)
            self.add_OFF_proper(off_proper)

        self._handle_excluded_propers()

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
            self.add_improper(PELE_improper)
            self.add_OFF_improper(off_improper)

    def add_atom(self, atom):
        """
        It adds an atom to the molecule's list of atoms.

        Parameters
        ----------
        atom : an peleffy.topology.Atom
            The Atom to add
        """
        self._atoms.append(atom)

    def add_bond(self, bond):
        """
        It adds a bond to the molecule's list of bonds.

        Parameters
        ----------
        bond : an peleffy.topology.Bond
            The Bond to add
        """
        self._bonds.append(bond)

    def add_angle(self, angle):
        """
        It adds an angle to the molecule's list of angles.

        Parameters
        ----------
        angle : an peleffy.topology.Angle
            The Angle to add
        """
        self._angles.append(angle)

    def add_proper(self, proper):
        """
        It adds a proper dihedral to the molecule's list of propers.

        Parameters
        ----------
        proper : an peleffy.topology.Proper
            The Proper to add
        """
        self._propers.append(proper)

    def add_OFF_proper(self, proper):
        """
        It adds a proper dihedral to the molecule's list of OFF propers.

        Parameters
        ----------
        proper : an peleffy.topology.OFFProper
            The OFFProper to add
        """
        self._OFF_propers.append(proper)

    def add_improper(self, improper):
        """
        It adds an improper dihedral to the molecule's list of impropers.

        Parameters
        ----------
        improper : an peleffy.topology.Improper
            The Improper to add
        """
        self._impropers.append(improper)

    def add_OFF_improper(self, improper):
        """
        It adds an improper dihedral to the molecule's list of OFF impropers.

        Parameters
        ----------
        improper : an peleffy.topology.OFFImproper
            The OFFImproper to add
        """
        self._OFF_impropers.append(improper)

    @property
    def molecule(self):
        """
        It returns the molecule that belongs to this topology element.

        Returns
        -------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object belonging to the topology
        """
        return self._molecule

    @property
    def parameters(self):
        """
        It returns the parameter wrapper that belongs to this topology
        element.

        Returns
        -------
        parameters : a BaseParameterWrapper
            The parameter wrapper belonging to the topology
        """
        return self._parameters

    @property
    def atoms(self):
        """
        The list of atoms in the topology.

        Returns
        -------
        atoms : list[peleffy.topology.molecule.Atom]
            The list of atoms of this Topology object.
        """
        return self._atoms

    @property
    def bonds(self):
        """
        The list of bonds in the topology.

        Returns
        -------
        bonds : list[peleffy.topology.Bond]
            The list of bonds of this Topology object.
        """
        return self._bonds

    @property
    def angles(self):
        """
        The list of angles in the topology.

        Returns
        -------
        angles : list[peleffy.topology.Angle]
            The list of angles of this Topology object.
        """
        return self._angles

    @property
    def propers(self):
        """
        The list of propers in the topology.

        Returns
        -------
        propers : list[peleffy.topology.Proper]
            The list of propers of this Topology object.
        """
        return self._propers

    @property
    def impropers(self):
        """
        The list of impropers in the topology.

        Returns
        -------
        impropers : list[peleffy.topology.Improper]
            The list of impropers of this Topology object.
        """
        return self._impropers
