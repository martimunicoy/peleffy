"""
This module handles all classes and functions related with the mapping
of two molecular topologies.
"""


class Mapper(object):
    """
    It defines the Mapper class.
    """

    _TIMEOUT = 150  # Timeout to find the MCS, in seconds

    def __init__(self, molecule1, molecule2, include_hydrogens=True):
        """
        Given two molecules, it finds the maximum common substructure
        (MCS) and maps their atoms.

        Parameters
        ----------
        molecule1: a peleffy.topology.Molecule
            The first molecule to map
        molecule2: a peleffy.topology.Molecule
            The second molecule to map
        include_hydrogens: bool
            Whether to include hydrogen atoms in the mapping or not.
            Default is True
        """

        # Check parameters
        import peleffy

        if (not isinstance(molecule1, peleffy.topology.Molecule)
                and not
                isinstance(molecule1, peleffy.topology.molecule.Molecule)):
            raise TypeError('Invalid input molecule 1')

        if (not isinstance(molecule2, peleffy.topology.Molecule)
                and not
                isinstance(molecule2, peleffy.topology.molecule.Molecule)):
            raise TypeError('Invalid input molecule 2')

        if molecule1.rdkit_molecule is None:
            raise ValueError('Molecule 1 has not been initialized')

        if molecule2.rdkit_molecule is None:
            raise ValueError('Molecule 2 has not been initialized')

        self._molecule1 = molecule1
        self._molecule2 = molecule2
        self._include_hydrogens = include_hydrogens

    def get_mcs(self):
        """
        It returns the Maximum Common Substructure (MCS) between
        both molecules.

        Parameters
        ----------
        mcs_mol : an RDKit.molecule object
            The resulting MCS molecule
        """
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        rdkit_toolkit = RDKitToolkitWrapper()

        mcs_mol = rdkit_toolkit.get_mcs(self.molecule1, self.molecule2,
                                        self._include_hydrogens,
                                        self._TIMEOUT)

        return mcs_mol

    def get_mapping(self):
        """
        It returns the mapping between both molecules.

        Returns
        -------
        mapping : list[tuple]
            The list of atom pairs between both molecules, represented
            with tuples
        """
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        rdkit_toolkit = RDKitToolkitWrapper()

        mcs_mol = self.get_mcs()

        mapping = rdkit_toolkit.get_atom_mapping(self.molecule1,
                                                 self.molecule2,
                                                 mcs_mol,
                                                 self._include_hydrogens)

        return mapping

    @property
    def molecule1(self):
        """
        It returns the first molecule to map.

        Returns
        -------
        molecule1 : a peleffy.topology.Molecule
            The first molecule to map
        """
        return self._molecule1

    @property
    def molecule2(self):
        """
        It returns the second molecule to map.

        Returns
        -------
        molecule2 : a peleffy.topology.Molecule
            The second molecule to map
        """
        return self._molecule2

    def to_png(self, output_png):
        """
        It generates a PNG image representing the resulting alchemical
        mapping.

        Parameters
        ----------
        output_png : str
            Path to the output PNG file to write
        """
        import os
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        extension = os.path.splitext(output_png)[1]

        if extension != ".png":
            raise ValueError("Invalid extension for a PNG file")


        rdkit_toolkit = RDKitToolkitWrapper()

        mcs_mol = rdkit_toolkit.get_mcs(self.molecule1, self.molecule2,
                                        self._include_hydrogens,
                                        self._TIMEOUT)

        image = rdkit_toolkit.draw_mapping(self.molecule1, self.molecule2,
                                           mcs_mol, self._include_hydrogens)

        image.save(output_png)

    def _ipython_display_(self):
        """
        It returns a representation of the mapping.

        Returns
        -------
        mapping_representation : a IPython display object
            Displayable RDKit molecules with mapping information
        """
        from IPython.display import display
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        rdkit_toolkit = RDKitToolkitWrapper()

        mcs_mol = rdkit_toolkit.get_mcs(self.molecule1, self.molecule2,
                                        self._include_hydrogens,
                                        self._TIMEOUT)

        image = rdkit_toolkit.draw_mapping(self.molecule1, self.molecule2,
                                           mcs_mol, self._include_hydrogens)

        return display(image)
