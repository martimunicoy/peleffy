"""
This module contains classes and methods related with alchemical modifications for molecular topologies.
"""


from abc import ABC
from dataclasses import dataclass


class Alchemizer(object):
    """
    It defines the Alchemizer class.
    """

    def __init__(self, topology1, topology2):
        """
        It initializes an Alchemizer object, which generates alchemical
        representations considering the topologies of two different
        molecules.

        Parameters
        ----------
        topology1 : a peleffy.topology.Topology object
            The molecular topology representation of molecule 1
        topology2 : a peleffy.topology.Topology object
            The molecular topology representation of molecule 2
        """

        # Check topologies
        import peleffy

        if (not isinstance(topology1, peleffy.topology.Topology)
                and not
                isinstance(topology1, peleffy.topology.topology.Topology)):
            raise TypeError('Invalid input topology 1')

        if (not isinstance(topology2, peleffy.topology.Topology)
                and not
                isinstance(topology2, peleffy.topology.topology.Topology)):
            raise TypeError('Invalid input topology 2')

        self._topology1 = topology1
        self._topology2 = topology2
        self._molecule1 = topology1.molecule
        self._molecule2 = topology2.molecule

        from peleffy.topology import Mapper

        # Map atoms from both molecules
        mapper = Mapper(self.molecule1, self.molecule2,
                        include_hydrogens=True)

        self._mapping = mapper.get_mapping()
        self._mcs_mol = mapper.get_mcs()

        # Join the two topologies
        self._joint_topology, self._non_native_atoms, \
            self._non_native_bonds, self._non_native_angles, \
            self._non_native_propers, self._non_native_impropers = \
            self._join_topologies()

        # Get exclusive topological elements in topology 1
        self._exclusive_atoms, self._exclusive_bonds, \
        self._exclusive_angles, self._exclusive_propers, \
        self._exclusive_impropers = self._get_exclusive_elements()

    @property
    def topology1(self):
        """
        It returns the first peleffy's Topology.

        Returns
        -------
        topology1 : a peleffy.topology.Topology
            The first peleffy's Topology object
        """
        return self._topology1

    @property
    def topology2(self):
        """
        It returns the second peleffy's Topology.

        Returns
        -------
        topology2 : a peleffy.topology.Topology
            The second peleffy's Topology object
        """
        return self._topology2

    @property
    def molecule1(self):
        """
        It returns the first molecule.

        Returns
        -------
        molecule1 : a peleffy.topology.Molecule
            The first Molecule object
        """
        return self._molecule1

    @property
    def molecule2(self):
        """
        It returns the second molecule.

        Returns
        -------
        molecule2 : a peleffy.topology.Molecule
            The second Molecule object
        """
        return self._molecule2

    @property
    def mapping(self):
        """
        It returns the mapping between both molecules.

        Returns
        -------
        mapping : list[tuple]
            The list of atom pairs between both molecules, represented
            with tuples
        """
        return self._mapping

    @property
    def mcs_mol(self):
        """
        It returns the Maximum Common Substructure (MCS) between
        both molecules.

        Parameters
        ----------
        mcs_mol : an RDKit.molecule object
            The resulting MCS molecule
        """
        return self._mcs_mol

    def get_alchemical_topology(self, fep_lambda=None, coul1_lambda=None,
                                coul2_lambda=None, vdw_lambda=None,
                                bonded_lambda=None):
        """
        Given a lambda, it returns an alchemical topology after
        combining both input topologies.

        Parameters
        ----------
        fep_lambda : float
            The value to define an FEP lambda. This lambda affects
            all the parameters. It needs to be contained between
            0 and 1. Default is None
        coul1_lambda : float
            The value to define a coulombic lambda for exclusive atoms
            of molecule 1. This lambda only affects coulombic parameters
            of exclusive atoms of molecule 1. It needs to be contained
            between 0 and 1. Default is None
        coul2_lambda : float
            The value to define a coulombic lambda for exclusive atoms
            of molecule 2. This lambda only affects coulombic parameters
            of exclusive atoms of molecule 2. It needs to be contained
            between 0 and 1. Default is None
        vdw_lambda : float
            The value to define a vdw lambda. This lambda only
            affects van der Waals parameters. It needs to be contained
            between 0 and 1. Default is None
        bonded_lambda : float
            The value to define a coulombic lambda. This lambda only
            affects bonded parameters. It needs to be contained
            between 0 and 1. Default is None

        Returns
        -------
        alchemical_topology : a peleffy.topology.Topology
            The resulting alchemical topology
        """
        # Define lambdas
        fep_lambda = FEPLambda(fep_lambda)
        coul1_lambda = Coulombic1Lambda(coul1_lambda)
        coul2_lambda = Coulombic2Lambda(coul2_lambda)
        vdw_lambda = VanDerWaalsLambda(vdw_lambda)
        bonded_lambda = BondedLambda(bonded_lambda)

        lambda_set = LambdaSet(fep_lambda, coul1_lambda, coul2_lambda,
                               vdw_lambda, bonded_lambda)

        alchemical_topology = self.topology_from_lambda_set(lambda_set)

        return alchemical_topology

    def topology_from_lambda_set(self, lambda_set):
        """
        Given a lambda, it returns an alchemical topology after
        combining both input topologies.

        Parameters
        ----------
        lambda_set : a peleffy.topology.alchemy.LambdaSet object
            The set of lambdas to use in the generation of the
            alchemical topology

        Returns
        -------
        alchemical_topology : a peleffy.topology.Topology
            The resulting alchemical topology
        """
        from copy import deepcopy

        alchemical_topology = deepcopy(self._joint_topology)

        for atom_idx, atom in enumerate(alchemical_topology.atoms):
            if atom_idx in self._exclusive_atoms:
                atom.apply_lambda(["sigma", "epsilon", "born_radius",
                                   "SASA_radius", "nonpolar_gamma",
                                   "nonpolar_alpha"],
                                  lambda_set.get_lambda_for_vdw(),
                                  reverse=False)
                atom.apply_lambda(["charge"],
                                  lambda_set.get_lambda_for_coulomb1(),
                                  reverse=False)

            if atom_idx in self._non_native_atoms:
                atom.apply_lambda(["sigma", "epsilon", "born_radius",
                                   "SASA_radius", "nonpolar_gamma",
                                   "nonpolar_alpha"],
                                  lambda_set.get_lambda_for_vdw(),
                                  reverse=True)
                atom.apply_lambda(["charge"],
                                  lambda_set.get_lambda_for_coulomb2(),
                                  reverse=True)

        for bond_idx, bond in enumerate(alchemical_topology.bonds):
            if bond_idx in self._exclusive_bonds:
                bond.apply_lambda(["spring_constant"],
                                  lambda_set.get_lambda_for_bonded(),
                                  reverse=False)

            if bond_idx in self._non_native_bonds:
                bond.apply_lambda(["spring_constant"],
                                  lambda_set.get_lambda_for_bonded(),
                                  reverse=True)

        for angle_idx, angle in enumerate(alchemical_topology.angles):
            if angle_idx in self._exclusive_angles:
                angle.apply_lambda(["spring_constant"],
                                   lambda_set.get_lambda_for_bonded(),
                                   reverse=False)

            if angle_idx in self._non_native_angles:
                angle.apply_lambda(["spring_constant"],
                                   lambda_set.get_lambda_for_bonded(),
                                   reverse=True)

        for proper_idx, proper in enumerate(alchemical_topology.propers):
            if proper_idx in self._exclusive_propers:
                proper.apply_lambda(["spring_constant"],
                                    lambda_set.get_lambda_for_bonded(),
                                    reverse=False)

            if proper_idx in self._non_native_propers:
                proper.apply_lambda(["spring_constant"],
                                    lambda_set.get_lambda_for_bonded(),
                                    reverse=True)

        for improper_idx, improper in enumerate(alchemical_topology.impropers):
            if improper_idx in self._exclusive_propers:
                improper.apply_lambda(["spring_constant"],
                                      lambda_set.get_lambda_for_bonded(),
                                      reverse=False)

            if improper_idx in self._non_native_propers:
                improper.apply_lambda(["spring_constant"],
                                      lambda_set.get_lambda_for_bonded(),
                                      reverse=True)

        return alchemical_topology

    def _join_topologies(self):
        """
        It joins the both topologies into a single one that contains
        all topological elements.

        Returns
        -------
        joint_topology : a peleffy.topology.Topology
            The resulting alchemical topology
        non_native_atoms : list[int]
            The list of atom indices that were added to topology 1
        non_native_bonds : list[int]
            The list of bond indices that were added to topology 1
        non_native_angles : list[int]
            The list of angle indices that were added to topology 1
        non_native_propers : list[int]
            The list of proper indices that were added to topology 1
        non_native_impropers : list[int]
            The list of improper indices that were added to topology 1
        """
        from copy import deepcopy
        from peleffy.topology import Topology
        from peleffy.utils.toolkits import RDKitToolkitWrapper
        from peleffy.utils import Logger

        # General initializers
        logger = Logger()
        rdkit_wrapper = RDKitToolkitWrapper()

        # First initialize the joint topology with topology 1
        joint_topology = deepcopy(self.topology1)

        # Initialize list of non native topological elements
        non_native_atoms = list()
        non_native_bonds = list()
        non_native_angles = list()
        non_native_propers = list()
        non_native_impropers = list()

        # Define mappers
        mol1_mapped_atoms = [atom_pair[0] for atom_pair in self.mapping]
        mol2_mapped_atoms = [atom_pair[1] for atom_pair in self.mapping]
        mol2_to_mol1_map = dict(zip(mol2_mapped_atoms, mol1_mapped_atoms))

        # Add atoms from topology 2 that are missing in topology 1
        atom_names = self.molecule1.get_pdb_atom_names()
        mol2_elements = rdkit_wrapper.get_elements(self.molecule2)
        mol2_to_alc_map = dict()
        for atom_idx, atom in enumerate(self.topology2.atoms):
            if atom_idx not in mol2_mapped_atoms:
                new_atom = deepcopy(atom)

                # Handle index of new atom
                new_index = len(joint_topology.atoms)
                new_atom.set_index(new_index)
                mol2_to_alc_map[atom_idx] = new_index

                # Handle name of new atom
                counter = 1
                new_name = '{:^4}'.format(str(mol2_elements[atom_idx]) +
                                          str(counter))
                while new_name in atom_names:
                    counter += 1
                    new_name = '{:^4}'.format(str(mol2_elements[atom_idx]) +
                                              str(counter))
                atom_names.append(new_name)
                new_atom.set_PDB_name(new_name)

                # Handle core allocation of new atom
                new_atom.set_as_branch()

                # Add new atom to the alchemical topology
                non_native_atoms.append(new_index)
                joint_topology.add_atom(new_atom)

            else:
                mol2_to_alc_map[atom_idx] = mol2_to_mol1_map[atom_idx]

        # Get an arbitrary absolute parent for molecule 2 (as long as it is in the MCS)
        absolute_parent = None
        for atom_idx in range(0, len(self.topology2.atoms)):
            if atom_idx in mol2_mapped_atoms:
                absolute_parent = atom_idx
                break
        else:
            logger.error(['Error: no atom in the MCS found in molecule ' +
                          f'{self.molecule1.name}'])

        # Get parent ids according to this new absolute parent of molecule 2
        mol2_parent_idxs = self.molecule2.graph.get_parents(absolute_parent)

        # Assign parents
        for atom_idx, atom in enumerate(self.topology2.atoms):
            if atom_idx not in mol2_mapped_atoms:
                alc_child_idx = mol2_to_alc_map[atom_idx]
                alc_parent_idx = mol2_to_alc_map[mol2_parent_idxs[atom_idx]]
                alc_parent_atom = joint_topology.atoms[alc_parent_idx]
                joint_topology.atoms[alc_child_idx].set_parent(
                    alc_parent_atom)

        # Add bonds
        for bond in self.topology2.bonds:
            atom1_idx = bond.atom1_idx
            atom2_idx = bond.atom2_idx

            if (atom1_idx not in mol2_mapped_atoms or
                    atom2_idx not in mol2_mapped_atoms):
                new_bond = deepcopy(bond)

                # Handle index of new atom
                new_index = len(joint_topology.bonds)
                new_bond.set_index(new_index)

                # Handle atom indices
                new_bond.set_atom1_idx(mol2_to_alc_map[atom1_idx])
                new_bond.set_atom2_idx(mol2_to_alc_map[atom2_idx])

                # Add new bond to the alchemical topology
                joint_topology.add_bond(new_bond)

        # Add angles
        for angle in self.topology2.angles:
            atom1_idx = angle.atom1_idx
            atom2_idx = angle.atom2_idx
            atom3_idx = angle.atom3_idx

            if (atom1_idx not in mol2_mapped_atoms or
                atom2_idx not in mol2_mapped_atoms or
                    atom3_idx not in mol2_mapped_atoms):
                new_angle = deepcopy(angle)

                # Handle index of new atom
                new_index = len(joint_topology.angles)
                new_angle.set_index(new_index)

                # Handle atom indices
                new_angle.set_atom1_idx(mol2_to_alc_map[atom1_idx])
                new_angle.set_atom2_idx(mol2_to_alc_map[atom2_idx])
                new_angle.set_atom3_idx(mol2_to_alc_map[atom3_idx])

                # Add new bond to the alchemical topology
                non_native_angles.append(new_index)
                joint_topology.add_angle(new_angle)

        # Add propers
        for proper in self.topology2.propers:
            atom1_idx = proper.atom1_idx
            atom2_idx = proper.atom2_idx
            atom3_idx = proper.atom3_idx
            atom4_idx = proper.atom1_idx

            if (atom1_idx not in mol2_mapped_atoms or
                atom2_idx not in mol2_mapped_atoms or
                atom3_idx not in mol2_mapped_atoms or
                    atom4_idx not in mol2_mapped_atoms):
                new_proper = deepcopy(proper)

                # Handle index of new atom
                new_index = len(joint_topology.propers)
                new_proper.set_index(new_index)

                # Handle atom indices
                new_proper.set_atom1_idx(mol2_to_alc_map[atom1_idx])
                new_proper.set_atom2_idx(mol2_to_alc_map[atom2_idx])
                new_proper.set_atom3_idx(mol2_to_alc_map[atom3_idx])
                new_proper.set_atom4_idx(mol2_to_alc_map[atom4_idx])

                # Add new bond to the alchemical topology
                non_native_propers.append(new_index)
                joint_topology.add_proper(new_proper)

        # Add impropers
        for improper in self.topology2.impropers:
            atom1_idx = improper.atom1_idx
            atom2_idx = improper.atom2_idx
            atom3_idx = improper.atom3_idx
            atom4_idx = improper.atom1_idx

            if (atom1_idx not in mol2_mapped_atoms or
                atom2_idx not in mol2_mapped_atoms or
                atom3_idx not in mol2_mapped_atoms or
                    atom4_idx not in mol2_mapped_atoms):
                new_improper = deepcopy(improper)

                # Handle index of new atom
                new_index = len(joint_topology.impropers)
                new_improper.set_index(new_index)

                # Handle atom indices
                new_improper.set_atom1_idx(mol2_to_alc_map[atom1_idx])
                new_improper.set_atom2_idx(mol2_to_alc_map[atom2_idx])
                new_improper.set_atom3_idx(mol2_to_alc_map[atom3_idx])
                new_improper.set_atom4_idx(mol2_to_alc_map[atom4_idx])

                # Add new bond to the alchemical topology
                non_native_impropers.append(new_index)
                joint_topology.add_improper(new_improper)

        return joint_topology, non_native_atoms, non_native_bonds, \
            non_native_angles, non_native_propers, non_native_impropers

    def _get_exclusive_elements(self):
        """
        It identifies those topological elements that are exclusive
        of topology 1. The condition to be exclusive is to belong to
        topology 1 but not to the MCS.

        Returns
        -------
        exclusive_atoms : list[int]
            The list of atom indices that are exclusive of topology 1
        exclusive_bonds : list[int]
            The list of bond indices that are exclusive of topology 1
        exclusive_angles : list[int]
            The list of angle indices that are exclusive of topology 1
        exclusive_propers : list[int]
            The list of proper indices that are exclusive of topology 1
        exclusive_impropers : list[int]
            The list of improper indices that are exclusive of topology 1
        """
        # Define mappers
        mol1_mapped_atoms = [atom_pair[0] for atom_pair in self.mapping]

        # Initialize list of exclusive topological elements
        exclusive_atoms = list()
        exclusive_bonds = list()
        exclusive_angles = list()
        exclusive_propers = list()
        exclusive_impropers = list()

        # Identify atoms
        for atom_idx, atom in enumerate(self.topology1.atoms):
            if atom_idx not in mol1_mapped_atoms:
                exclusive_atoms.append(atom_idx)

        # Identify bonds
        for bond_idx, bond in enumerate(self.topology1.bonds):
            atom1_idx = bond.atom1_idx
            atom2_idx = bond.atom2_idx

            if (atom1_idx not in mol1_mapped_atoms or
                    atom2_idx not in mol1_mapped_atoms):
                exclusive_bonds.append(bond_idx)

        # Identify angles
        for angle_idx, angle in enumerate(self.topology1.angles):
            atom1_idx = angle.atom1_idx
            atom2_idx = angle.atom2_idx
            atom3_idx = angle.atom3_idx

            if (atom1_idx not in mol1_mapped_atoms or
                atom2_idx not in mol1_mapped_atoms or
                    atom3_idx not in mol1_mapped_atoms):
                exclusive_angles.append(angle_idx)

       # Identify propers
        for proper_idx, proper in enumerate(self.topology1.propers):
            atom1_idx = proper.atom1_idx
            atom2_idx = proper.atom2_idx
            atom3_idx = proper.atom3_idx
            atom4_idx = proper.atom3_idx

            if (atom1_idx not in mol1_mapped_atoms or
                atom2_idx not in mol1_mapped_atoms or
                atom3_idx not in mol1_mapped_atoms or
                    atom4_idx not in mol1_mapped_atoms):
                exclusive_propers.append(proper_idx)

       # Identify impropers
        for improper_idx, improper in enumerate(self.topology1.impropers):
            atom1_idx = improper.atom1_idx
            atom2_idx = improper.atom2_idx
            atom3_idx = improper.atom3_idx
            atom4_idx = improper.atom3_idx

            if (atom1_idx not in mol1_mapped_atoms or
                atom2_idx not in mol1_mapped_atoms or
                atom3_idx not in mol1_mapped_atoms or
                    atom4_idx not in mol1_mapped_atoms):
                exclusive_impropers.append(improper_idx)

        return exclusive_atoms, exclusive_bonds, exclusive_angles, \
               exclusive_propers, exclusive_impropers


class Lambda(ABC):
    """
    It defines the Lambda class.
    """

    _TYPE = ""

    def __init__(self, value=None):
        """
        It initializes a Lambda object.

        Parameters
        ----------
        value : float
            The value of this Lambda object. It needs to be
            contained between 0 and 1. Default is None
        """
        if value is not None:
            try:
                value = float(value)
            except ValueError:
                raise ValueError("Invalid value for a lambda: " +
                                 f"\'{value}\'")
            if (value > 1) or (value < 0):
                raise ValueError("Invalid value for a lambda: " +
                                 f"\'{value}\'. " +
                                 "It has to be between 0 and 1")

        self._value = value
        pass

    @property
    def value(self):
        """
        It returns the value of this Lambda object.

        Returns
        -------
        value : float or None
            The value of this Lambda object. It can be None if the
            value for this Lambda object has not been set
        """
        return self._value

    @property
    def type(self):
        """
        It returns the type of this Lambda object.

        Returns
        -------
        type : str
            The type of this Lambda object
        """
        return self._TYPE

    @property
    def is_set(self):
        """
        It answers whether the value of this Lambda object has
        been set or not.

        Returns
        -------
        is_set : bool
            It is true only if the value of Lambda object has been
            set
        """
        return self.value is not None


class FEPLambda(Lambda):
    """
    It defines the FEPLambda class. It affects all parameters.
    """
    _TYPE = "fep"


class Coulombic1Lambda(Lambda):
    """
    It defines the CoulombicLambda1 class. It affects only coulombic
    parameters involving exclusive atoms of molecule 1.
    """
    _TYPE = "coulombic1"


class Coulombic2Lambda(Lambda):
    """
    It defines the CoulombicLambda2 class. It affects only coulombic
    parameters involving exclusive atoms of molecule 2.
    """
    _TYPE = "coulombic2"


class VanDerWaalsLambda(Lambda):
    """
    It defines the VanDerWaalsLambda class. It affects only van der Waals
    parameters.
    """
    _TYPE = "vdw"


class BondedLambda(Lambda):
    """
    It defines the BondedLambda class. It affects only bonded parameters.
    """
    _TYPE = "bonded"


class LambdaSet(object):
    """
    It defines the LambdaSet class.
    """

    def __init__(self, fep_lambda, coul1_lambda, coul2_lambda,
                 vdw_lambda, bonded_lambda):
        """
        It initializes a LambdaSet object which stores all the different
        types of lambda.

        Parameters
        ----------
        fep_lambda : a peleffy.topology.alchemy.FEPLambda object
            The fep lambda
        coul1_lambda : a peleffy.topology.alchemy.CoulombicLambda object
            The coulombic lambda for exclusive atoms of molecule 1
        coul2_lambda : a peleffy.topology.alchemy.CoulombicLambda object
            The coulombic lambda for exclusive atoms of molecule 2
        vdw_lambda : a peleffy.topology.alchemy.VanDerWaalsLambda object
            The van der Waals lambda
        bonded_lambda : a peleffy.topology.alchemy.BondedLambda object
            The bonded lambda
        """
        import peleffy

        # Check parameters
        if not isinstance(fep_lambda,
                          peleffy.topology.alchemy.FEPLambda):
            raise TypeError('Invalid fep_lambda supplied to LambdaSet')
        if not isinstance(coul1_lambda,
                          peleffy.topology.alchemy.Coulombic1Lambda):
            raise TypeError('Invalid coul1_lambda supplied to LambdaSet')
        if not isinstance(coul2_lambda,
                          peleffy.topology.alchemy.Coulombic2Lambda):
            raise TypeError('Invalid coul2_lambda supplied to LambdaSet')
        if not isinstance(vdw_lambda,
                          peleffy.topology.alchemy.VanDerWaalsLambda):
            raise TypeError('Invalid vdw_lambda supplied to LambdaSet')
        if not isinstance(bonded_lambda,
                          peleffy.topology.alchemy.BondedLambda):
            raise TypeError('Invalid bonded_lambda supplied to LambdaSet')

        self._fep_lambda = fep_lambda
        self._coul1_lambda = coul1_lambda
        self._coul2_lambda = coul2_lambda
        self._vdw_lambda = vdw_lambda
        self._bonded_lambda = bonded_lambda

    @property
    def fep_lambda(self):
        """
        It returns the fep_lambda value.

        Returns
        -------
        fep_lambda : float
            The value of the fep_lambda
        """
        return self._fep_lambda

    @property
    def coul1_lambda(self):
        """
        It returns the coul1_lambda value.

        Returns
        -------
        coul1_lambda : float
            The value of the coul1_lambda
        """
        return self._coul1_lambda

    @property
    def coul2_lambda(self):
        """
        It returns the coul2_lambda value.

        Returns
        -------
        coul2_lambda : float
            The value of the coul2_lambda
        """
        return self._coul2_lambda

    @property
    def vdw_lambda(self):
        """
        It returns the vdw_lambda value.

        Returns
        -------
        vdw_lambda : float
            The value of the vdw_lambda
        """
        return self._vdw_lambda

    @property
    def bonded_lambda(self):
        """
        It returns the bonded_lambda value.

        Returns
        -------
        bonded_lambda : float
            The value of the bonded_lambda
        """
        return self._bonded_lambda

    def get_lambda_for_vdw(self):
        """
        It returns the lambda to be applied on van der Waals parameters.

        Returns
        -------
        lambda_value : float
            The lambda value to be applied on van der Waals parameters
        """
        if self.vdw_lambda.is_set:
            lambda_value = self.vdw_lambda.value

        elif self.fep_lambda.is_set:
            lambda_value = self.fep_lambda.value

        else:
            lambda_value = 0.0

        return lambda_value

    def get_lambda_for_coulomb1(self):
        """
        It returns the lambda to be applied on Coulomb parameters of
        exclusive atoms of molecule 1.

        Returns
        -------
        lambda_value : float
            The lambda value to be applied on Coulomb parameters of
            exclusive atoms of molecule 1
        """
        if self.coul1_lambda.is_set:
            lambda_value = self.coul1_lambda.value

        elif self.fep_lambda.is_set:
            lambda_value = self.fep_lambda.value

        else:
            lambda_value = 0.0

        return lambda_value

    def get_lambda_for_coulomb2(self):
        """
        It returns the lambda to be applied on Coulomb parameters of
        exclusive atoms of molecule 2.

        Returns
        -------
        lambda_value : float
            The lambda value to be applied on Coulomb parameters of
            exclusive atoms of molecule 2
        """
        if self.coul2_lambda.is_set:
            lambda_value = self.coul2_lambda.value

        elif self.fep_lambda.is_set:
            lambda_value = self.fep_lambda.value

        else:
            lambda_value = 0.0

        return lambda_value

    def get_lambda_for_bonded(self):
        """
        It returns the lambda to be applied on bonded parameters.

        Returns
        -------
        lambda_value : float
            The lambda value to be applied on bonded parameters
        """
        if self.bonded_lambda.is_set:
            lambda_value = self.bonded_lambda.value

        elif self.fep_lambda.is_set:
            lambda_value = self.fep_lambda.value

        else:
            lambda_value = 0.0

        return lambda_value
