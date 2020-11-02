"""
This module contains classes and functions involved in the manipulation of
PELE's solvent templates.
"""

from simtk import unit

from peleffy.utils import get_data_file_path, warning_on_one_line
from peleffy.utils import Logger


class _SolventWrapper(object):
    """
    A wrapper for any solvent-like class.
    """
    _ff_file = None
    _name = None

    def __init__(self, molecule):
        """
        Initializes a SolventWrapper object.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        self._molecule = molecule
        self._radii = dict.fromkeys([tuple((idx, ))
                                     for idx in range(0, len(molecule.atoms))],
                                    unit.Quantity())
        self._scales = dict.fromkeys([tuple((idx, ))
                                      for idx in range(0, len(molecule.atoms))],
                                     unit.Quantity())
        self._solvent_dielectric = float(0)
        self._solute_dielectric = float(0)
        self._surface_area_penalty = float(0)
        self._solvent_radius = float(0)
        self._initialize_from_molecule()

    def _initialize_from_molecule(self):
        """
        Initializes a SolventWrapper object using an peleffy's Molecule.
        """
        logger = Logger()
        logger.info(' - Loading solvent parameters')

        from peleffy.utils.toolkits import OpenForceFieldToolkitWrapper

        off_toolkit = OpenForceFieldToolkitWrapper()
        GBSA_handler = off_toolkit.get_parameter_handler_from_forcefield(
            'GBSA', self._ff_file)

        self._solvent_dielectric = GBSA_handler.solvent_dielectric
        self._solute_dielectric = GBSA_handler.solute_dielectric
        self._surface_area_penalty = GBSA_handler.surface_area_penalty
        self._solvent_radius = GBSA_handler.solvent_radius

        from peleffy.forcefield import OpenForceField

        forcefield = OpenForceField(self._ff_file)
        parameters = forcefield.parameterize(self.molecule)

        self._radii = parameters['GBSA_radii']
        self._scales = parameters['GBSA_scales']

    def to_dict(self):
        """
        Returns this SolventWrapper object as a dictionary.

        Returns
        -------
        data : dict
            A dictionary containing the data of this SolventWrapper object
        """
        data = dict()
        data['SolventParameters'] = dict()
        data['SolventParameters']['Name'] = self.name
        data['SolventParameters']['General'] = dict()
        data['SolventParameters']['General']['solvent_dielectric'] = \
            round(self.solvent_dielectric, 5)
        data['SolventParameters']['General']['solute_dielectric'] = \
            round(self.solute_dielectric, 5)
        data['SolventParameters']['General']['solvent_radius'] = \
            round(self.solvent_radius.value_in_unit(unit.angstrom), 5)
        data['SolventParameters']['General']['surface_area_penalty'] = \
            round(self.surface_area_penalty.value_in_unit(
                unit.kilocalorie / (unit.angstrom**2 * unit.mole)), 8)
        data['SolventParameters'][self.molecule.tag] = dict()

        atom_names = self.molecule.get_pdb_atom_names()

        for atom, name in zip(self.molecule.rdkit_molecule.GetAtoms(),
                              atom_names):
            name = name.replace(' ', '_')
            index = atom.GetIdx()
            data['SolventParameters'][self.molecule.tag][name] = \
                {'radius': round(self.radii[tuple((index, ))].value_in_unit(
                                 unit.angstrom), 5),
                 'scale': round(self.scales[tuple((index, ))], 5)}

        return data

    def to_json_file(self, path):
        """
        Writes this SolventWrapper object to a json file.

        Parameters
        ----------
        path : str
            Path to save the json file to
        """
        import json
        with open(path, 'w') as file:
            json.dump(self.to_dict(), file, indent=4)

    @property
    def name(self):
        """
        The name of the solvent.

        Returns
        -------
        name : str
            The name of this solvent object.
        """
        return self._name

    @property
    def molecule(self):
        """
        The peleffy's Molecule to parameterize.

        Returns
        -------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        return self._molecule

    @property
    def radii(self):
        """
        The dict of radii of the parameterized molecule.

        Returns
        -------
        radii : dict[atom indexes: simtk.unit.Quantity]
            The radius assigned to each atom of the molecule
        """
        return self._radii

    @property
    def scales(self):
        """
        The dict of scales of the parameterized molecule.

        Returns
        -------
        scales : dict[atom indexes: simtk.unit.Quantity]
            The scale assigned to each atom of the molecule
        """
        return self._scales

    @property
    def solvent_dielectric(self):
        """
        The solvent dielectric value of this solvent object.

        Returns
        -------
        solvent_dielectric : float
            The solvent dielectric value
        """
        return self._solvent_dielectric

    @property
    def solute_dielectric(self):
        """
        The solute dielectric value of this solvent object.

        Returns
        -------
        solute_dielectric : float
            The solute dielectric value
        """
        return self._solute_dielectric

    @property
    def surface_area_penalty(self):
        """
        The surface area penalty value of this solvent object.

        Returns
        -------
        surface_area_penalty : float
            The surface area penalty value
        """
        return self._surface_area_penalty

    @property
    def solvent_radius(self):
        """
        The solvent radius value of this solvent object.

        Returns
        -------
        solvent_radius : float
            The solvent radius value
        """
        return self._solvent_radius


class OBC1(_SolventWrapper):
    """
    Implementation of the OBC1 solvent.
    """

    _ff_file = get_data_file_path('forcefields/GBSA_OBC1-1.0.offxml')
    _name = 'OBC1'

    def __init__(self, molecule):
        """
        Initializes an OBC1 object.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        # Not implemented in PELE
        import warnings
        warnings.formatwarning = warning_on_one_line
        warnings.warn("OBC1 is not implemented in PELE", Warning)

        super().__init__(molecule)

    def _initialize_from_molecule(self):
        """
        Initializes the OBC1 solvent using an peleffy's Molecule.
        """
        super()._initialize_from_molecule()


class OBC2(_SolventWrapper):
    """
    Implementation of the OBC2 solvent.
    """

    _ff_file = get_data_file_path('forcefields/GBSA_OBC2-1.0.offxml')
    _name = 'OBC2'

    def __init__(self, molecule):
        """
        Initializes an OBC2 object.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file

        Examples
        --------

        Generate the solvent parameters of a molecule

        >>> from peleffy.topology import Molecule
        >>> from peleffy.solvent import OBC2

        >>> molecule = Molecule('molecule.pdb')
        >>> solvent = OBC2(molecule)
        >>> solvent.to_json_file('molecule_solv.json')

        """
        super().__init__(molecule)

    def _initialize_from_molecule(self):
        """
        Initializes the OBC2 solvent using an peleffy's Molecule.
        """
        super()._initialize_from_molecule()
