import numpy as np
from simtk import unit

from offPELE.utils import get_data_file_path, warning_on_one_line
from offPELE.utils.toolkits import OpenForceFieldToolkitWrapper


class _SolventWrapper(object):
    _ff_file = None
    _name = None

    def __init__(self, molecule):
        self._molecule = molecule
        self._radii = np.zeros(len(self.molecule.atoms))
        self._scales = np.zeros(len(self.molecule.atoms))
        self._solvent_dielectric = float(0)
        self._solute_dielectric = float(0)
        self._surface_area_penalty = float(0)
        self._solvent_radius = float(0)
        self._initialize_from_molecule()

    def _initialize_from_molecule(self):
        off_toolkit = OpenForceFieldToolkitWrapper()

        GBSA_handler = off_toolkit.get_parameter_handler_from_forcefield(
            'GBSA', self._ff_file)

        self._solvent_dielectric = GBSA_handler.solvent_dielectric
        self._solute_dielectric = GBSA_handler.solute_dielectric
        self._surface_area_penalty = GBSA_handler.surface_area_penalty
        self._solvent_radius = GBSA_handler.solvent_radius

        parameters = off_toolkit.get_parameters_from_forcefield(
            self._ff_file, self.molecule)

        self._radii = parameters.get_GBSA_radii()
        self._scales = parameters.get_GBSA_scales()

    def to_dict(self):
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
        data['SolventParameters'][self.molecule.name] = dict()
        for atom in self.molecule.rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            name = pdb_info.GetName().replace(' ', '_')
            index = atom.GetIdx()
            data['SolventParameters'][self.molecule.name][name] = \
                {'radius': round(self.radii[tuple((index, ))].value_in_unit(
                                 unit.angstrom), 5),
                 'scale': round(self.scales[tuple((index, ))], 5)}

        return data

    def to_json_file(self, path):
        import json
        with open(path, 'w') as file:
            json.dump(self.to_dict(), file, indent=4)

    @property
    def name(self):
        return self._name

    @property
    def molecule(self):
        return self._molecule

    @property
    def radii(self):
        return self._radii

    @property
    def scales(self):
        return self._scales

    @property
    def solvent_dielectric(self):
        return self._solvent_dielectric

    @property
    def solute_dielectric(self):
        return self._solute_dielectric

    @property
    def surface_area_penalty(self):
        return self._surface_area_penalty

    @property
    def solvent_radius(self):
        return self._solvent_radius


class OBC1(_SolventWrapper):
    _ff_file = get_data_file_path('forcefields/GBSA_OBC1-1.0.offxml')
    _name = 'OBC1'

    def __init__(self, molecule):
        # Not implemented in PELE
        import warnings
        warnings.formatwarning = warning_on_one_line
        warnings.warn("OBC1 is not implemented in PELE", Warning)

        super().__init__(molecule)

    def _initialize_from_molecule(self):
        super()._initialize_from_molecule()


class OBC2(_SolventWrapper):
    _ff_file = get_data_file_path('forcefields/GBSA_OBC1-1.0.offxml')
    _name = 'OBC2'

    def __init__(self, molecule):
        super().__init__(molecule)

    def _initialize_from_molecule(self):
        super()._initialize_from_molecule()
