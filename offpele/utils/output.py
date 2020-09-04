"""
This module contains a set of classes and methods designed to handle
output workflows.
"""


class OutputPathHandler(object):
    """
    It handles the output paths of offpele parameter files.
    """
    IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OFF/Parsley/HeteroAtoms/'
    ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'
    SOLVENT_TEMPLATE_PATH = 'DataLocal/OBC/'

    def __init__(self, molecule, as_DataLocal=False):
        """
        It initializes an OutputPathHandler object.

        Parameters
        ----------
        molecule : An offpele.topology.Molecule
            A Molecule object to be written as an Impact file
        as_DataLocal : bool
            Whether to save output files following PELE's DataLocal
            hierarchy or not
        """
        self._molecule = molecule
        self._as_DataLocal = as_DataLocal

    def set_DataLocal_behaviour(self, as_DataLocal):
        """
        It sets the DataLocal behaviour of this OutputPathHandler object.

        Parameters
        ----------
        as_DataLocal : bool
            Whether to save output files following PELE's DataLocal
            hierarchy or not
        """
        self._as_DataLocal = as_DataLocal

    @property
    def as_DataLocal(self):
        """
        The DataLocal configuration of this OutputPathHandler object.

        Returns
        -------
        as_DataLocal : bool
            Whether to save output files following PELE's DataLocal
            hierarchy or not
        """
        return self._as_DataLocal

    @property
    def molecule(self):
        """
        The offpele's Molecule.

        Returns
        -------
        molecule : an offpele.topology.Molecule
            The offpele's Molecule object
        """
        return self._molecule
