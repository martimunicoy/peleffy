"""
This module contains a set of classes and methods designed to handle
output workflows.
"""


import os

from peleffy.utils import create_path


class OutputPathHandler(object):
    """
    It handles the output paths of peleffy parameter files.
    """

    OFF_IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OFF/Parsley/HeteroAtoms/'
    OPLS_IMPACT_TEMPLATE_PATH = 'DataLocal/Templates/OPLS2005/HeteroAtoms/'
    ROTAMER_LIBRARY_PATH = 'DataLocal/LigandRotamerLibs/'
    SOLVENT_TEMPLATE_PATH = 'DataLocal/OBC/'
    FILE_TYPES = ['impact template', 'rotamer library', 'solvent template']

    def __init__(self, molecule, as_datalocal=False, output_path=None):
        """
        It initializes an OutputPathHandler object.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        as_datalocal : bool
            Whether to save output files following PELE's DataLocal
            hierarchy or not
        output_path : str
            The output path inside which all files will be saved. Default
            is the current directory
        """
        self._molecule = molecule
        self._as_datalocal = as_datalocal
        if output_path is None:
            self._output_path = os.path.curdir
        else:
            self._output_path = output_path

    def set_DataLocal_behaviour(self, as_datalocal):
        """
        It sets the DataLocal behaviour of this OutputPathHandler object.

        Parameters
        ----------
        as_datalocal : bool
            Whether to save output files following PELE's DataLocal
            hierarchy or not
        """
        self._as_datalocal = as_datalocal

    def get_path(self, file_type, create_missing_folders=True):
        """
        It returns the path to a file, according to the file type.

        Parameters
        ----------
        file_type : str
            The file type whose path is requested. One of
            ['impact template', 'rotamer library', 'solvent template']
        create_missing_folders : bool
            Whether to create missing folders or not. Default is True

        Returns
        -------
        file_path : str
            The path for the requested file type
        """
        if file_type.lower() not in self.FILE_TYPES:
            raise ValueError('Invalid file type: {}'.format(file_type)
                             + ', it must be one of', self.FILE_TYPES)

        if file_type.lower() == 'impact template':
            return self.get_impact_template_path(create_missing_folders)

        if file_type.lower() == 'rotamer library':
            return self.get_rotamer_library_path(create_missing_folders)

        if file_type.lower() == 'solvent template':
            return self.get_solvent_template_path(create_missing_folders)

    def get_impact_template_path(self, create_missing_folders=True):
        """
        It returns the path for an Impact template file.

        Parameters
        ----------
        create_missing_folders : bool
            Whether to create missing folders or not. Default is True

        Returns
        -------
        file_path : str
            The path for an Impact template file
        """
        file_name = self._molecule.tag.lower() + 'z'

        if self.as_datalocal:
            forcefield = self.molecule.forcefield
            if forcefield is not None:
                if forcefield.type == 'OpenFF':
                    path = os.path.join(self.output_path,
                                        self.OFF_IMPACT_TEMPLATE_PATH)
                elif forcefield.type == 'OpenFF + OPLS2005':
                    if forcefield._nonbonding == 'openff':
                        path = os.path.join(self.output_path,
                                            self.OFF_IMPACT_TEMPLATE_PATH)
                    else:
                        path = os.path.join(self.output_path,
                                            self.OPLS_IMPACT_TEMPLATE_PATH)
                else:
                    path = os.path.join(self.output_path,
                                        self.OPLS_IMPACT_TEMPLATE_PATH)
        else:
            path = self.output_path

        if create_missing_folders:
            create_path(path)

        return os.path.join(path, file_name)

    def get_rotamer_library_path(self, create_missing_folders=True):
        """
        It returns the path for a rotamer library file.

        Parameters
        ----------
        create_missing_folders : bool
            Whether to create missing folders or not. Default is True

        Returns
        -------
        file_path : str
            The path for a rotamer library file
        """
        file_name = self._molecule.tag.upper() + '.rot.assign'

        if self.as_datalocal:
            path = os.path.join(self.output_path, self.ROTAMER_LIBRARY_PATH)
        else:
            path = self.output_path

        if create_missing_folders:
            create_path(path)

        return os.path.join(path, file_name)

    def get_solvent_template_path(self, create_missing_folders=True):
        """
        It returns the path for a solvent template file.

        Parameters
        ----------
        create_missing_folders : bool
            Whether to create missing folders or not. Default is True

        Returns
        -------
        file_path : str
            The path for a solvent template file
        """
        file_name = 'ligandParams.txt'

        if self.as_datalocal:
            path = os.path.join(self.output_path, self.SOLVENT_TEMPLATE_PATH)
        else:
            path = self.output_path

        if create_missing_folders:
            create_path(path)

        return os.path.join(path, file_name)

    @property
    def as_datalocal(self):
        """
        The DataLocal configuration of this OutputPathHandler object.

        Returns
        -------
        as_datalocal : bool
            Whether to save output files following PELE's DataLocal
            hierarchy or not
        """
        return self._as_datalocal

    @property
    def molecule(self):
        """
        The peleffy's Molecule.

        Returns
        -------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        return self._molecule

    @property
    def output_path(self):
        """
        The output path inside which all files will be saved.

        Returns
        output_path : str
            The output path
        """
        return self._output_path
