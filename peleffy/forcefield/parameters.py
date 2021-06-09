"""
This module provides classes to define and construct force field
parameters.
"""


__all__ = ["OpenForceFieldParameterWrapper",
           "OPLS2005ParameterWrapper",
           "OpenFFOPLS2005ParameterWrapper"]


from collections import defaultdict
from simtk import unit

from peleffy.utils import get_data_file_path
from peleffy.utils import Logger


class BaseParameterWrapper(dict):
    """
    It defines the base class of a parameter wrapper that inherits from
    dict.
    """

    _name = ''
    _fixed_keys = ('atom_names', 'atom_types', 'charges', 'sigmas',
                   'epsilons', 'SGB_radii', 'vdW_radii', 'gammas',
                   'alphas', 'GBSA_radii', 'GBSA_scales')
    _unfixed_keys = ('bonds', 'angles', 'propers', 'impropers')
    _attribute_restrictions = {'bonds': ('atom1_idx', 'atom2_idx',
                                         'spring_constant', 'eq_dist'),
                               'angles': ('atom1_idx', 'atom2_idx',
                                          'atom3_idx', 'spring_constant',
                                          'eq_angle'),
                               'propers': ('atom1_idx', 'atom2_idx',
                                           'atom3_idx', 'atom4_idx',
                                           'periodicity', 'phase',
                                           'k', 'idivf'),
                               'impropers': ('atom1_idx', 'atom2_idx',
                                             'atom3_idx', 'atom4_idx',
                                             'periodicity', 'phase',
                                             'k', 'idivf')}
    _keys = _fixed_keys + _unfixed_keys

    def __init__(self, parameters_dict=dict(), forcefield_name=None):
        """
        It initializes a parameter wrapper object.

        Parameters
        ----------
        parameters_dict : dict
             A dictionary keyed by parameter type that contains the
             parameters to store. Default is an empty dict
        forcefield_name : str
            The name of the force field the parameters belong to
        """
        for key in self._keys:
            self[key] = list()

        for key, value in parameters_dict.items():
            self[key] = value

        if forcefield_name is None:
            self._forcefield_name = self._name
        else:
            self._forcefield_name = forcefield_name

    def __eq__(self, other):
        """
        It sets the equality operator for the BaseParameterWrapper class.

        Parameters
        ----------
        other : a BaseParameterWrapper object
            The other BaseParameterWrapper object to compare with the
            current one
        """
        if not isinstance(other, BaseParameterWrapper):
            return False

        return super().__eq__(other) and \
            self.forcefield_name == other.forcefield_name

    def __ne__(self, other):
        """
        It sets the inequality operator for the BaseParameterWrapper class.

        Parameters
        ----------
        other : a BaseParameterWrapper object
            The other BaseParameterWrapper object to compare with the
            current one
        """
        return not self.__eq__(other)

    def __setitem__(self, key, val):
        """
        It sets an item in the dictionary, only if the key is an
        expected key and fulfills the size and attribute restrictions.

        Parameters
        ----------
        key : str
            The parameter type
        val : list[float, int or simtk.unit.Quantity]
            The list of corresponding values
        """
        if key not in self._keys:
            raise KeyError('An unexpected parameter type cannot be set '
                           + 'as a dictionary key. Unexpected type: '
                           + '{}'.format(key))

        okay = True
        if key in self._fixed_keys:
            for local_key in self.keys():
                if local_key not in self._fixed_keys:
                    continue
                local_size = len(self[local_key])
                val_size = len(val)
                if local_size != 0 and val_size != 0:
                    if len(val) != local_size:
                        okay = False

        if not okay:
            raise ValueError('Supplied values for key \'{}\' '.format(key)
                             + 'are not of the same size as values that '
                             + 'have already been loaded. Input length: '
                             + '{} Expected length: '.format(len(val))
                             + '{}'.format(local_size))

        okay = True
        if key in self._attribute_restrictions:
            attributes = self._attribute_restrictions[key]
            if isinstance(val, dict):
                if len(attributes) != len(val.keys()):
                    okay = False
                for attribute in attributes:
                    if attribute not in val.keys():
                        okay = False
            elif isinstance(val, list):
                for val_element in val:
                    if len(attributes) != len(val_element.keys()):
                        okay = False
                    for attribute in attributes:
                        if attribute not in val_element.keys():
                            okay = False
            else:
                raise ValueError('Unexpected input value with key '
                                 + '\'{}\' '.format(key)
                                 + 'and value \'{}\''.format(val))
        if not okay:
            raise ValueError('Supplied values for key \'{}\' '.format(key)
                             + 'do not contain the expected attributes. '
                             + 'They are: '
                             + '\'{}\''.format(', '.join(attributes)))

        dict.__setitem__(self, key, val)

    def __str__(self):
        """
        It returns the string representation of this parameter wrapper.

        Returns
        -------
        string_representation : str
            The string representation
        """
        return self.to_string()

    def add_parameters(self, label, parameters):
        """
        It adds a list of parameters of the same type to the collection.

        Parameters
        ----------
        label : str
            The label to describe the parameter type
        parameters : list
            The list of parameters to include to the collection
        """
        self[label] = parameters

    def to_string(self):
        """
        It returns a string representation of the parameters stored in
        the wrapper.

        Returns
        -------
        string_repr : str
            The string representation
        """

        from peleffy.utils import convert_all_quantities_to_string
        import pprint

        string_repr = pprint.pformat(convert_all_quantities_to_string(self))

        return string_repr

    def to_json(self, output_path):
        """
        It saves this parameter wrapper as a json file.

        Parameters
        ----------
        output_path : str
            The path to save the output json file
        """

        import json
        from peleffy.utils import convert_all_quantities_to_string

        with open(output_path, 'w') as f:
            json.dump(convert_all_quantities_to_string(self), f,
                      indent=4, sort_keys=True)

    def from_json(self, input_path):
        """
        It loads a json file containing the parameters into the parameter
        wrapper.

        Parameters
        ----------
        input_path : str
            The path to the json file to load the parameters from

        Returns
        -------
        self : a BaseParameterWrapper object
            The resulting parameters wrapper
        """
        import json
        def correct_type(label, value):
            """
            It converts the parameters loaded from the JSON file into the
            expected data type in the parameter wrapper.

            Parameters
            ----------
            label : str
                Label for the parameters
            value : list
                Set of parameters

            Returns
            -------
            correct_value : list
                Set of paramaeters with the expected data type
            """
            import simtk.unit
            from peleffy.utils import string_to_quantity

            # Dictionary relating parameters key with the expected data type
            dict_units = {
                'alphas': float, 'gammas': float,
                'charges':  simtk.unit.quantity.Quantity,
                'sigmas' : 'list_Quantity', 'epsilons' : 'list_Quantity',
                'SGB_radii' : 'list_Quantity','vdW_radii' : 'list_Quantity',
                'angles': dict, 'bonds': dict, 'impropers': dict,
                'propers': dict, 'GBSA_radii': 'list_Quantity',
                'GBSA_scales': 'list_Quantity',
                'atom_names': str, 'atom_types' : str
                }

            # Skip data type transformation if the list is empty or None values
            if value == [None,] * len(value):
                return value
            if value == []:
                return None

            # Perform data type transformation, if needed
            else:
                if dict_units[label] is str:
                    correct_value = value

                if dict_units[label] is float:
                    correct_value  = [float(v) if not v is None else v
                    for v in value]

                if dict_units[label] is simtk.unit.quantity.Quantity:
                    correct_value = string_to_quantity(value)
                    return correct_value

                if dict_units[label] == 'list_Quantity':
                    correct_value = \
                            [unit.Quantity(value=string_to_quantity(v)._value,
                                           unit=string_to_quantity(v).unit)
                            for v in value]

                if dict_units[label] is dict:
                    correct_value = []
                    for element in value:
                        correct_dict = dict()
                        for k,v in element.items():
                            if 'idx' in k:
                                correct_dict[k] = int(v)
                            else:
                                try:
                                    correct_dict[k] = unit.Quantity(
                                            value=string_to_quantity(v)._value,
                                            unit=string_to_quantity(v).unit)
                                except:
                                    correct_dict[k] = v
                        correct_value.append(correct_dict)
                return correct_value

        # Load the dict from the JSON file
        with open(input_path) as f:
            data = json.load(f)

        # Correct the data type format and fetch the BaseParameterWrapper
        for key, value in data.items():
            value_correct = correct_type(key, value)
            if not value_correct is None:
                self.add_parameters(key, value_correct)
        return self

    @property
    def atom_iterator(self):
        """
        It returns an iterator that goes through all atom parameters.
        The order of parameters that is followed is: 'atom_names',
        'atom_types', 'sigmas', 'epsilons', 'charges', 'SGB_radii',
        'vdW_radii', 'gammas', 'alphas'.

        Returns
        -------
        iterator : a zip object
            It is an iterator that returns all atom parameters following
            their order
        """
        assert len(self['atom_names']) == len(self['atom_types']) \
            and len(self['atom_names']) == len(self['sigmas']) \
            and len(self['atom_names']) == len(self['epsilons']) \
            and len(self['atom_names']) == len(self['charges']) \
            and len(self['atom_names']) == len(self['SGB_radii']) \
            and len(self['atom_names']) == len(self['vdW_radii']) \
            and len(self['atom_names']) == len(self['gammas']) \
            and len(self['atom_names']) == len(self['alphas']), \
            'Size of atom parameter lists should match'

        return zip(self['atom_names'], self['atom_types'],
                   self['sigmas'], self['epsilons'],
                   self['charges'], self['SGB_radii'],
                   self['vdW_radii'], self['gammas'],
                   self['alphas'])

    def is_empty(self):
        """
        If the are no parameters stored in this wrapper, it returns True.

        Returns
        -------
        answer : bool
            Whether this parameter wrapper is empty or not
        """
        return all(len(self[key]) == 0 for key in self)

    @property
    def name(self):
        """
        It returns the name of the current parameter wrapper.

        Returns
        -------
        name : str
            The name of the parameter wrapper
        """
        return self._name

    @staticmethod
    def from_impact_template(molecule, impact_template_path):
        """
        It returns a parameter wrapper out of an impact template.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object
        impact_template_path : str
            The path to the impact template from where the parameters
            will be fetched

        Returns
        -------
        params : a BaseParameterWrapper object
            The resulting parameters wrapper

        Examples
        --------

        Obtain a parameter wrapper from an Impact template.

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')
        >>> impact_template_path = 'molz'

        >>> from peleffy.forcefield.parameters import \
                OpenForceFieldParametersWrapper

        >>> wrapper_off = OpenForceFieldParametersWrapper()
        >>> parameters = wrapper_off.from_impact_template(
                molecule, impact_template_path)

        """

        def index(atom_idx):
            """
            Atom's index.

            Parameters
            ----------
            index : int
                The Atom's index for a WritableWrapper

            Returns
            -------
            index : int
                The Atom's index for a BaseParameterWrapper object
            """
            return abs(int(atom_idx)) - 1

        def zero_to_none(values):
            """
            It converts a vector of zeros to None.

            Parameters
            ----------
            values : list
                List of parameters.

            Returns
            -------
            values : list
                It is set to None if all the parameters are zeros.
            """
            n_atoms = len(values)
            if type(values[0]) == float:
                if all(value == float(0.0) for value in values):
                    values = [None, ] * n_atoms
            if type(values[0]) == unit.quantity.Quantity:
                if all(value == unit.Quantity(0, unit.angstroms) for value
                       in values):
                    values = values = [None, ] * n_atoms
            return values

        def get_phase(info):
            """
            It computes the phase of a Dihedral given its prefactor. PELE
            prefactors can be  equal to 1 or -1. If the phase constant is
            different than 0 or 180, an extra column in the Impact template
            specifies its value.

            Parameters
            ----------
            info : list
                List of parameters for this Dihedral

            Returns
            -------
            phase : simtk.unit.Quantity
                The phase constant of the Dihedral
            """
            if len(info) > 7:
                return unit.Quantity(
                    value=float(info[7]),
                    unit=unit.degree)
            else:
                prefactor = float(info[5])
                if prefactor == 1:
                    return unit.Quantity(
                        value=0.00,
                        unit=unit.degree)
                if prefactor == -1:
                    return unit.Quantity(
                        value=180.00,
                        unit=unit.degree)

        def reindex_atom_idx(list_params, dict_index):
            """
            It sorts and reindexes atoms of a parameters list according to a
            dictionary.

            Parameters
            ----------
            list_params : list
                List of parameters to reindex
            dict_index : dict
                Dictionary that defines the reindexer between old and new index

            Returns
            -------
            list_params : list
                Updated list of parameters with the correct indexes
            """
            for i, atom_info in enumerate(list_params):
                for key, value in atom_info.items():
                    if '_idx' in key:
                        list_params[i][key] = dict_index.get(value)
            return list_params

        # Parse Impact template file
        tag, nbon, bond, thet, phi, iphi = ([] for i in range(6))
        with open(impact_template_path, 'r') as fd:
            type_info = 'tag'
            for line in fd.readlines():
                if len(line.strip()) == 0:  # Skip empty lines
                    continue
                if not line.startswith('*'):
                    if 'NBON' in line:
                        type_info = 'nbon'
                    elif 'BOND' in line:
                        type_info = 'bond'
                    elif 'THET' in line:
                        type_info = 'thet'
                    elif 'PHI' in line and 'IPHI' not in line:
                        type_info = 'phi'
                    elif 'IPHI' in line:
                        type_info = 'iphi'
                    elif 'END' in line:
                        break
                    else:
                        if type_info == 'tag':
                            tag.append(line)
                        if type_info == 'nbon':
                            nbon.append(line)
                        if type_info == 'bond':
                            bond.append(line)
                        if type_info == 'thet':
                            thet.append(line)
                        if type_info == 'phi':
                            phi.append(line)
                        if type_info == 'iphi':
                            iphi.append(line)

        # Atom names and atom types
        atom_names_list, atom_types_list = ([] for i in range(2))
        for line in tag[1:]:
            info = line.split()
            atom_names_list.append(info[4])
            atom_types_list.append(info[3])

        # Charges, epsilons, sigmas, SGB radius, vdW radius, gammas and alphas
        charges_list, epsilons_list, sigmas_list, SGB_list, vdW_list, \
            gammas_list, alphas_list = ([] for i in range(7))
        for line in nbon:
            info = line.split()
            charges_list.append(unit.Quantity(
                value=float(info[3]),
                unit=unit.elementary_charge))
            epsilons_list.append(unit.Quantity(
                value=float(info[2]),
                unit=unit.kilocalorie / unit.mole))
            sigmas_list.append(unit.Quantity(
                value=float(info[1]),
                unit=unit.angstrom))
            SGB_list.append(unit.Quantity(
                value=float(info[4]),
                unit=unit.angstrom))
            vdW_list.append(unit.Quantity(
                value=float(info[5]),
                unit=unit.angstrom))
            gammas_list.append(float(info[6]))
            alphas_list.append(float(info[7]))

        # Bonds
        bonds_list = []
        for line in bond:
            info = line.split()
            case = {'atom1_idx': index(info[0]),
                    'atom2_idx': index(info[1]),
                    'eq_dist': unit.Quantity(
                        value=float(info[3]),
                        unit=unit.angstrom),
                    'spring_constant': unit.Quantity(
                        value=float(info[2]),
                        unit=unit.kilocalorie
                        / (unit.angstrom ** 2 * unit.mole))}
            bonds_list.append(case)

        # Angles
        angles_list = []
        for line in thet:
            info = line.split()
            case = {'atom1_idx': index(info[0]),
                    'atom2_idx': index(info[1]),
                    'atom3_idx': index(info[2]),
                    'eq_angle': unit.Quantity(
                        value=float(info[4]),
                        unit=unit.degrees),
                    'spring_constant': unit.Quantity(
                        value=float(info[3]),
                        unit=unit.kilocalorie
                        / (unit.radian ** 2 * unit.mole))}
            angles_list.append(case)

        # Impropers
        impropers_list = []
        for line in iphi:
            info = line.split()
            case = {'atom1_idx': index(info[0]),
                    'atom2_idx': index(info[1]),
                    'atom3_idx': index(info[2]),
                    'atom4_idx': index(info[3]),
                    'idivf': abs(float(info[5])),
                    'k': unit.Quantity(
                        value=float(info[4]),
                        unit=unit.kilocalorie / unit.mole),
                    'periodicity': int(float(info[6])),
                    'phase': get_phase(info)}
            impropers_list.append(case)

        # Propers
        propers_list = []
        for line in phi:
            info = line.split()
            case = {'atom1_idx': index(info[0]),
                    'atom2_idx': index(info[1]),
                    'atom3_idx': index(info[2]),
                    'atom4_idx': index(info[3]),
                    'idivf': abs(float(info[5])),
                    'k': unit.Quantity(
                        value=float(info[4]),
                        unit=unit.kilocalorie / unit.mole),
                    'periodicity': int(float(info[6])),
                    'phase': get_phase(info)}
            propers_list.append(case)

        # Initialize the BaseParameterWrapper object
        params = BaseParameterWrapper()

        # PDB atom names from Molecule object
        ordered_pdb_atom_names = [pdb_name.replace(' ', '_') for pdb_name in
                                  molecule.get_pdb_atom_names()]

        # Check up that the Molecule representation and the Impact template
        # represent the same chemical entity
        if not set(ordered_pdb_atom_names) == set(atom_names_list):
            raise ValueError(
                "The Impact template file {} ".format(impact_template_path) +
                "does not represent the same chemical entity as the molecule.")

        # Molecule object and Impact template have the atoms in the same order
        if ordered_pdb_atom_names == atom_names_list:
            # Assign parameters from Impact template to the BaseParameterWrapper
            # object
            params['angles'] = angles_list
            params['atom_names'] = atom_names_list
            params['atom_types'] = atom_types_list
            params['bonds'] = bonds_list
            params['charges'] = charges_list
            params['epsilons'] = epsilons_list
            params['impropers'] = impropers_list
            params['propers'] = propers_list
            params['sigmas'] = sigmas_list
            params['vdW_radii'] = vdW_list
            params['SGB_radii'] = zero_to_none(SGB_list)
            params['gammas'] = zero_to_none(gammas_list)
            params['alphas'] = zero_to_none(alphas_list)

        # Molecule object and Impact template have the atoms in different order
        if ordered_pdb_atom_names != atom_names_list:
            # Reindexer dictionary to sort atoms as Molecule object
            index_order = [atom_names_list.index(atom_name) for atom_name in
                           ordered_pdb_atom_names]
            dict_reindex = dict(zip(index_order, range(len(index_order))))

            # Assign parameters from Impact template to the BaseParameterWrapper
            # object reordering atoms according to Molecule object atoms
            params['atom_names'] = [atom_names_list[i] for i in index_order]
            params['atom_types'] = [atom_types_list[i] for i in index_order]
            params['charges'] = [charges_list[i] for i in index_order]
            params['epsilons'] = [epsilons_list[i] for i in index_order]
            params['sigmas'] = [sigmas_list[i] for i in index_order]
            params['vdW_radii'] = [vdW_list[i] for i in index_order]
            params['SGB_radii'] = [zero_to_none(SGB_list)[i] for i in
                                   index_order]
            params['gammas'] = [zero_to_none(gammas_list)[i] for i in
                                index_order]
            params['alphas'] = [zero_to_none(alphas_list)[i] for i in
                                index_order]
            params['impropers'] = reindex_atom_idx(impropers_list,
                                                   dict_reindex)
            params['propers'] = reindex_atom_idx(propers_list,
                                                 dict_reindex)
            params['angles'] = reindex_atom_idx(angles_list,
                                                dict_reindex)
            params['bonds'] = reindex_atom_idx(bonds_list,
                                               dict_reindex)
        return params

    @property
    def forcefield_name(self):
        return self._forcefield_name


class OpenForceFieldParameterWrapper(BaseParameterWrapper):
    """
    It defines a parameters wrapper for an OpenFF force field.
    """

    _name = 'OpenFF'

    @staticmethod
    def from_label_molecules(molecule, parameters, openff_version):
        """
        It parses the parameters coming from the label_molecules()
        method from OpenFF toolkit.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        parameters_dict : dict
             A dictionary keyed by parameter type that contains the
             parameters to store
        openff_version : str
            The version of the OpenFF force field the parameters belong to

        Returns
        -------
        parameters : an OpenForceFieldParameterWrapper object
            The resulting parameters wrapper
        """

        I_6R_OF_2 = 0.8908987181403393  # The inverse of the sixth root of 2

        def build_dict(parameters, attribute_name):
            """
            It builds the dictionary of parameters of a specific force type.

            Parameters
            ----------
            parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters of a specific force type grouped by tuples
                with the atom ids that the parameters belong to
            attribute_name : str
                The name of the attribute that is requested

            Returns
            -------
            value_by_index : dict[tuple, parameter_value]
                The parameter values that were requested grouped by the atom
                ids they belong to (arranged as a tuple)
            """
            if parameters:
                value_by_index = dict()
                for index, parameter in parameters.items():
                    value_by_index[index] = getattr(parameter, attribute_name)

                return value_by_index

        def build_dynamic_dicts(parameters, attr_core_name):
            """
            It builds the dynamically the dictionaries of parameters of a
            specific force type.

            It works with the same idea as _build_dict(), however it can
            handle multiple dictionary definitions were consecutive
            parameters of the same type are found in the force type's
            parameters dictionary. It works, for example, with the multiple
            proper and improper definitions found in the OpenForceField API.
            More information in the <ProperTorsions> and <ImproperTorsions>
            sections at:
            https://open-forcefield-toolkit.readthedocs.io/en/0.7.0/smirnoff.html

            Parameters
            ----------
            parameters : dict[tuple, openforcefield.typing.engines.smirnoff.parameters.ParameterHandler]
                The parameters of a specific force type grouped by tuples
                with the atom ids that the parameters belong to
            attribute_name : str
                The name of the attribute that is requested

            Returns
            -------
            value_by_index : dict[tuple, parameter_value]
                The parameter values that were requested grouped by the atom
                ids they belong to (arranged as a tuple)
            """
            if parameters:
                parameters_by_index = defaultdict(dict)
                all_attr_ids_found = set()
                for index, parameter in parameters.items():
                    counter = int(1)
                    attr_name = attr_core_name + str(counter)
                    while(attr_name in parameter.to_dict()):
                        all_attr_ids_found.add(counter)
                        attr_value = getattr(parameter, attr_name)
                        parameters_by_index[index][counter] = attr_value
                        counter += 1
                        attr_name = attr_core_name + str(counter)

                output_list = list()
                for attr_id in sorted(all_attr_ids_found):
                    value_by_index = dict()
                    for index in parameters.keys():
                        if attr_id in parameters_by_index[index]:
                            value_by_index[index] = \
                                parameters_by_index[index][attr_id]
                        else:
                            value_by_index[index] = None

                    output_list.append(value_by_index)

                return output_list

        def sigma_from_rmin_half(rmin_half):
            """
            It converts rmin_half values to sigmas according to:
            http://ambermd.org/Questions/vdwequation.pdf
            """
            sigma = I_6R_OF_2 * 2 * rmin_half

            return sigma

        # This will be the dictionary that will contain the parsed data
        params = defaultdict(list)

        pdb_atom_names = molecule.get_pdb_atom_names()
        n_atoms = len(pdb_atom_names)

        # PELE needs underscores instead of whitespaces
        params['atom_names'] = [name.replace(' ', '_') for name
                                in pdb_atom_names]
        params['atom_types'] = ['OFFT', ] * n_atoms

        # TODO: Safety check of the parameters, if any atom has no parameters
        # it should complain here

        # van der Waals parameters
        if 'vdW' in parameters:
            vdW_parameters = parameters['vdW']
            sigma_by_idx = build_dict(vdW_parameters, 'sigma')
            for idx in range(0, n_atoms):
                params['sigmas'].append(sigma_by_idx[(idx, )])

            if all([sigma is None for sigma in params['sigmas']]):
                params['sigmas'] = list()
                rmin_half_by_idx = build_dict(vdW_parameters, 'rmin_half')
                for idx in range(0, n_atoms):
                    params['sigmas'].append(sigma_from_rmin_half(
                        rmin_half_by_idx[(idx, )]))

            epsilon_by_idx = build_dict(vdW_parameters, 'epsilon')
            for idx in range(0, n_atoms):
                params['epsilons'].append(epsilon_by_idx[(idx, )])

            params['SGB_radii'] = [None, ] * n_atoms
            params['vdW_radii'] = [sigma / 2.0 for sigma in params['sigmas']]
            params['gammas'] = [None, ] * n_atoms
            params['alphas'] = [None, ] * n_atoms

        # bond parameters
        if 'Bonds' in parameters:
            bond_parameters = parameters['Bonds']
            length_by_idx = build_dict(bond_parameters, 'length')
            k_by_idx = build_dict(bond_parameters, 'k')
            for idxs in length_by_idx.keys():
                idx1, idx2 = idxs
                params['bonds'].append(
                    {'atom1_idx': idx1,
                     'atom2_idx': idx2,
                     # PELE works with half of the OFF's spring
                     'spring_constant': k_by_idx[idxs] / 2.0,
                     'eq_dist': length_by_idx[idxs]
                     })

        # angle parameters
        if 'Angles' in parameters:
            angle_parameters = parameters['Angles']
            angle_by_idx = build_dict(angle_parameters, 'angle')
            k_by_idx = build_dict(angle_parameters, 'k')
            for idxs in angle_by_idx.keys():
                idx1, idx2, idx3 = idxs
                params['angles'].append(
                    {'atom1_idx': idx1,
                     'atom2_idx': idx2,
                     'atom3_idx': idx3,
                     # PELE works with half of the OFF's spring
                     'spring_constant': k_by_idx[idxs] / 2.0,
                     'eq_angle': angle_by_idx[idxs]
                     })

        # proper parameters
        if 'ProperTorsions' in parameters:
            proper_parameters = parameters['ProperTorsions']
            periodicities = build_dynamic_dicts(proper_parameters,
                                                'periodicity')
            phases = build_dynamic_dicts(proper_parameters, 'phase')
            ks = build_dynamic_dicts(proper_parameters, 'k')
            idivfs = build_dynamic_dicts(proper_parameters, 'idivf')
            if (periodicities is not None and phases is not None
                    and ks is not None and idivfs is not None):
                for periodicity_by_idx, phase_by_idx, k_by_idx, \
                        idivf_by_idx in zip(periodicities, phases, ks,
                                            idivfs):

                    assert periodicity_by_idx.keys() == phase_by_idx.keys() \
                        and periodicity_by_idx.keys() == k_by_idx.keys() \
                        and periodicity_by_idx.keys() == idivf_by_idx.keys(), \
                        'Unconsistent torsional parameter indexes. '\
                        + 'Keys should match.'

                    for idxs in periodicity_by_idx.keys():
                        idx1, idx2, idx3, idx4 = idxs

                        periodicity = periodicity_by_idx[idxs]
                        phase = phase_by_idx[idxs]
                        k = k_by_idx[idxs]
                        idivf = idivf_by_idx[idxs]

                        if (periodicity is not None and phase is not None
                                and k is not None and idivf is not None):
                            params['propers'].append(
                                {'atom1_idx': idx1,
                                 'atom2_idx': idx2,
                                 'atom3_idx': idx3,
                                 'atom4_idx': idx4,
                                 'periodicity': periodicity,
                                 'phase': phase,
                                 'k': k,
                                 'idivf': idivf
                                 })

        # improper parameters
        if 'ImproperTorsions' in parameters:
            improper_parameters = parameters['ImproperTorsions']
            periodicities = build_dynamic_dicts(improper_parameters,
                                                'periodicity')
            phases = build_dynamic_dicts(improper_parameters, 'phase')
            ks = build_dynamic_dicts(improper_parameters, 'k')
            idivfs = build_dynamic_dicts(improper_parameters, 'idivf')
            if (periodicities is not None and phases is not None
                    and ks is not None):
                # idivf is a optional parameter in OpenForceField
                if len(idivfs) == 0:
                    for period_by_index in periodicities:
                        idivfs.append(dict(zip(
                            period_by_index.keys(),
                            [1, ] * len(period_by_index.keys()))))

                for periodicity_by_idx, phase_by_idx, k_by_idx, \
                        idivf_by_idx in zip(periodicities, phases, ks,
                                            idivfs):

                    assert periodicity_by_idx.keys() == phase_by_idx.keys() \
                        and periodicity_by_idx.keys() == k_by_idx.keys() \
                        and periodicity_by_idx.keys() == idivf_by_idx.keys(), \
                        'Unconsistent torsional parameter indexes. '\
                        + 'Keys should match.'

                    for idxs in periodicity_by_idx.keys():
                        idx1, idx2, idx3, idx4 = idxs

                        periodicity = periodicity_by_idx[idxs]
                        phase = phase_by_idx[idxs]
                        k = k_by_idx[idxs]
                        idivf = idivf_by_idx[idxs]

                        if (periodicity is not None and phase is not None
                                and k is not None and idivf is not None):
                            params['impropers'].append(
                                {'atom1_idx': idx1,
                                 'atom2_idx': idx2,
                                 'atom3_idx': idx3,
                                 'atom4_idx': idx4,
                                 'periodicity': periodicity,
                                 'phase': phase,
                                 'k': k,
                                 'idivf': idivf
                                 })

        # GBSA parameters
        if 'GBSA' in parameters:
            gbsa_parameters = parameters['GBSA']
            params['GBSA_radii'] = build_dict(gbsa_parameters, 'radius')
            params['GBSA_scales'] = build_dict(gbsa_parameters, 'scale')

        return OpenForceFieldParameterWrapper(params, openff_version)


class OPLS2005ParameterWrapper(BaseParameterWrapper):
    """
    It defines a parameters wrapper for OPLS2005 force field.
    """

    _name = 'OPLS2005'

    @staticmethod
    def from_ffld_output(molecule, ffld_output):
        """
        It parses the parameters coming from the Schrodinger's
        ffld_server output file.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object
        ffld_output : str
            The path to the ffld_server's output file

        Returns
        -------
        parameters : an OpenForceFieldParameterWrapper object
            The resulting parameters wrapper
        """

        from simtk import unit
        from peleffy.utils.toolkits import RDKitToolkitWrapper

        params = defaultdict(list)

        # Assign atom names according to the input PDB file (if any)
        pdb_atom_names = molecule.get_pdb_atom_names()

        # PELE needs underscores instead of whitespaces
        params['atom_names'] = [name.replace(' ', '_') for name
                                in pdb_atom_names]

        section = 'out'
        name_to_index = dict()  # To pair atom names and indexes

        lines_iterator = iter(ffld_output.split('\n'))

        for line in lines_iterator:
            if line.startswith('OPLSAA FORCE FIELD TYPE ASSIGNED'):
                section = 'atoms'

                # Skip the next 3 lines
                next(lines_iterator)
                next(lines_iterator)
                next(lines_iterator)

            elif line.startswith(' Stretch'):
                section = 'bonds'

            elif line.startswith(' Bending'):
                section = 'angles'

            elif line.startswith(' proper Torsion'):
                section = 'propers'

            elif line.startswith(' improper Torsion'):
                section = 'impropers'

            elif line == '':
                continue

            elif section == 'atoms':
                if line.startswith('-'):
                    continue

                fields = line.split()
                assert len(fields) > 7, 'Unexpected number of fields ' \
                    + 'found at line {}'.format(line)

                name_to_index[line[0:8]] = len(name_to_index)

                params['atom_types'].append(fields[3])
                params['charges'].append(
                    unit.Quantity(float(fields[4]),
                                  unit.elementary_charge))
                params['sigmas'].append(
                    unit.Quantity(float(fields[5]),
                                  unit.angstrom))
                params['epsilons'].append(
                    unit.Quantity(float(fields[6]),
                                  unit.kilocalorie / unit.mole))

            elif section == 'bonds':
                fields = line.split()
                assert len(fields) > 4, 'Unexpected number of fields ' \
                    + 'found at line {}'.format(line)

                params['bonds'].append(
                    {'atom1_idx': name_to_index[line[0:8]],
                     'atom2_idx': name_to_index[line[8:16]],
                     'spring_constant': unit.Quantity(
                        float(fields[2]), unit.kilocalorie
                        / (unit.angstrom ** 2 * unit.mole)),
                     'eq_dist': unit.Quantity(float(fields[3]),
                                              unit.angstrom)
                     })

            elif section == 'angles':
                fields = line.split()
                assert len(fields) > 5, 'Unexpected number of fields ' \
                    + 'found at line {}'.format(line)

                params['angles'].append(
                    {'atom1_idx': name_to_index[line[0:8]],
                     'atom2_idx': name_to_index[line[8:16]],
                     'atom3_idx': name_to_index[line[16:24]],
                     'spring_constant': unit.Quantity(
                        float(fields[3]), unit.kilocalorie
                        / (unit.radian ** 2 * unit.mole)),
                     'eq_angle': unit.Quantity(float(fields[4]),
                                               unit.degrees)
                     })

            elif section == 'propers':
                fields = line.split()
                assert len(fields) > 9, 'Unexpected number of fields ' \
                    + 'found at line {}'.format(line)

                atom1_idx = name_to_index[line[0:8]]
                atom2_idx = name_to_index[line[8:16]]
                atom3_idx = name_to_index[line[16:24]]
                atom4_idx = name_to_index[line[24:32]]

                for k, periodicity, phase in zip(
                        fields[4:8], [1, 2, 3, 4],
                        [unit.Quantity(0.0, unit.degree),
                         unit.Quantity(180.0, unit.degree),
                         unit.Quantity(0.0, unit.degree),
                         unit.Quantity(180.0, unit.degree)]):
                    k = float(k)
                    if k != 0.0:
                        k = unit.Quantity(k, unit.kilocalorie / unit.mole)
                        params['propers'].append(
                            {'atom1_idx': atom1_idx,
                             'atom2_idx': atom2_idx,
                             'atom3_idx': atom3_idx,
                             'atom4_idx': atom4_idx,
                             'periodicity': periodicity,
                             'phase': phase,
                             'k': k / 2.0,  # PELE works with half of Schrodinger's force constant
                             'idivf': 1.0
                             })

                # In case all four k's are zero, we still need to include
                # the proper torsion to be used by PELE in the 1-4
                # interactions
                if all([float(k) == 0.0 for k in fields[4:8]]):
                    params['propers'].append(
                        {'atom1_idx': atom1_idx,
                         'atom2_idx': atom2_idx,
                         'atom3_idx': atom3_idx,
                         'atom4_idx': atom4_idx,
                         'periodicity': 1.0,
                         'phase': unit.Quantity(0.0, unit.degree),
                         'k': unit.Quantity(0.0,
                                            unit.kilocalorie / unit.mole),
                         'idivf': 1.0
                         })

            elif section == 'impropers':
                fields = line.split()
                assert len(fields) > 5, 'Unexpected number of fields ' \
                    + 'found at line {}'.format(line)

                k = unit.Quantity(float(fields[4]),
                                  unit.kilocalorie / unit.mole)

                params['impropers'].append(
                    {'atom1_idx': name_to_index[line[0:8]],
                     'atom2_idx': name_to_index[line[8:16]],
                     'atom3_idx': name_to_index[line[16:24]],
                     'atom4_idx': name_to_index[line[24:32]],
                     'periodicity': 2,
                     'phase': unit.Quantity(180.0, unit.degree),
                     'k': k / 2.0,  # PELE works with half of Schrodinger's force constant
                     'idivf': 1.0
                     })

        opls_parameters_wrapper = OPLS2005ParameterWrapper(params)
        OPLS2005ParameterWrapper._add_SGBNP_solvent_parameters(
            opls_parameters_wrapper)

        # Employ RDKit to extract atom degree and parent type lists
        wrapper = RDKitToolkitWrapper()
        atom_names = wrapper.get_atom_names(molecule)
        degree_by_name = dict(zip(atom_names,
                                  wrapper.get_atom_degrees(molecule)))
        parent_by_name = dict(zip(atom_names,
                                  wrapper.get_hydrogen_parents(molecule)))
        element_by_name = dict(zip(atom_names,
                                   wrapper.get_elements(molecule)))

        OPLS2005ParameterWrapper._add_GBSA_solvent_parameters(
            opls_parameters_wrapper, degree_by_name,
            parent_by_name, element_by_name)

        return opls_parameters_wrapper

    @staticmethod
    def _find_similar_atom_types(atom_type, tried):
        """
        It tries to find a similar atom type, skipping the ones that
        have already been tried. It uses the definitions from the
        similarity.param file.

        Parameters
        ----------
        atom_type : str
            The atom type from which similar atom types will be searched
        tried : list[str]
            The list of atom types that have already been tried and
            will be skipped

        Returns
        -------
        new_atom_type : str
            The most similar atom type that has been found, if any.
            Otherwise, it returns None
        """

        new_atom_type = None
        best_similarity = 0
        similarity_path = get_data_file_path(
            'parameters/similarity.param')

        with open(similarity_path) as f:
            for line in f:
                fields = line.split()
                assert len(fields) > 2, 'Unexpected number of fields ' \
                    + 'at line {}'.format(line)

                atom_type1, atom_type2, similarity = fields[0:3]
                if (atom_type == atom_type1
                        and float(similarity) > best_similarity
                        and atom_type2 not in tried):
                    best_similarity = float(similarity)
                    new_atom_type = atom_type2
                elif (atom_type == atom_type2
                        and float(similarity) > best_similarity
                        and atom_type1 not in tried):
                    best_similarity = float(similarity)
                    new_atom_type = atom_type1

        return new_atom_type

    @staticmethod
    def _add_SGBNP_solvent_parameters(OPLS_params):
        """
        It adds the SGBNP solvent parameters (used in the SGBNP solvent
        implemented in the OPLS2005 of PELE) to the OPLS parameters
        collection.

        Parameters
        ----------
        OPLS_params : an OPLS2005ParameterWrapper object
            The set of lists of parameters grouped by parameter type.
            Thus, the dictionary has the following keys: atom_names,
            atom_types, charges, sigmas, and epsilons. The following
            solvent parameters will be added to the collection: SGB_radii,
            vdW_radii, gammas, alphas
        """
        from simtk import unit

        solvent_data = dict()
        parameters_path = get_data_file_path(
            'parameters/f14_sgbnp.param')

        with open(parameters_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.split()
                assert len(fields) > 7, 'Unexpected line with less ' \
                    'than 8 fields at {}'.format(line)

                atom_type = fields[1]

                solvent_data[atom_type] = {
                    'SGB_radii': unit.Quantity(float(fields[4]),
                                               unit.angstrom),
                    'vdW_radii': unit.Quantity(float(fields[5]),
                                               unit.angstrom),
                    'gammas': float(fields[6]),
                    'alphas': float(fields[7])}

        parameters_to_add = defaultdict(list)
        tried = list()

        for atom_type in OPLS_params['atom_types']:
            parameters_found = False
            while(not parameters_found):
                if atom_type in solvent_data:
                    for label, value in solvent_data[atom_type].items():
                        parameters_to_add[label].append(value)
                    parameters_found = True

                else:
                    new_atom_type = \
                        OPLS2005ParameterWrapper._find_similar_atom_types(
                            atom_type, tried)
                    if new_atom_type is None:
                        atom_type = 'DF'  # Set it to default
                    else:
                        tried.append(new_atom_type)
                        atom_type = new_atom_type

        for label, params in parameters_to_add.items():
            OPLS_params.add_parameters(label, params)

    @staticmethod
    def _add_GBSA_solvent_parameters(OPLS_params, degree_by_name,
                                     parent_by_name,
                                     element_by_name):
        """
        It adds the GBSA solvent parameters (used in the OBC solvent
        implemented in the OPLS2005 of PELE) to the OPLS parameters
        collection.

        Parameters
        ----------
        OPLS_params : an OPLS2005ParameterWrapper object
            The set of lists of parameters grouped by parameter type.
            Thus, the dictionary has the following keys: atom_names,
            atom_types, charges, sigmas, and epsilons. The following
            solvent parameters will be added to the collection: SGB_radii,
            vdW_radii, gammas, alphas
        degree_by_name : dict
            The dictionary containing the number of bonds for
            each atom_name
        parent_by_name : dict
            The dictionary containing the element of the parent
            for hydrogen atoms, keyed by atom name of the child
        element_by_name : dict
            The dictionary containing the atom elements keyed by
            atom name
        """
        import re
        import json

        def _check_bonds(scale, atom_type, degree, parent):
            """
            It checks the number of bonds and the parent atom for
            terminal H. Besides, in some especified cases it updates
            the scale factor.

            Parameters
            ----------
            scale : str
                The scale factor of the current atom
            atom_type : str
                The atom type of the current atom
            degree : int
                The number of bonds of the current atom
            parent : str
                The element belonging to the parent atom, in case that the
                current atom is a hydrogen atom
            """

            if atom_type == 'H' and parent == 'O':
                scale = float(1.05)
            if atom_type == 'H' and parent == 'N':
                scale = float(1.15)
            if atom_type == 'C' and degree == 3:
                scale = float(1.875)
            if atom_type == 'C' and degree == 2:
                scale = float(1.825)
            if atom_type == 'N' and degree == 4:
                scale = float(1.625)
            if atom_type == 'N' and degree == 1:
                scale = float(1.60)
            if atom_type == 'O' and degree == 1:
                scale = float(1.48)

            return scale

        def _find_GBSA_parameters_according_to(atom_name, atom_type,
                                               degree, element, parent):
            """
            It computes the HTC radii and the Overlap factor for Heteroatoms.
            The parameters have been extracted from Tinker Molecular package. If
            one parameter has not been defined in the OBC templates it puts the
            default parameter and raises a warning.

            Parameters
            ----------
            atom_name : str
                Atom name
            atom_type : str
                Atom type
            degree : str
                Number of bonds
            element : str
                The element the atom belongs to
            parent : str
                The element of the atom's parent (if the atom is a hydrogen)

            Returns
            -------
            radius : str
                HCT radii
            scale : str
                Overlap factor
            """

            PARAMS_PATH = get_data_file_path('parameters/OBCparam.json')

            # Get rid of terminal white spaces
            atom_name = atom_name.strip()
            atom_type = atom_type.strip()

            # Load the dictioaries with the OBC parameters
            with open(PARAMS_PATH) as fd:
                params_by_type, scale_by_element, \
                    radius_by_element = json.load(fd)

            # Assign scale factor and HCT radius using the atom type
            if atom_type in params_by_type:
                scale, radius = params_by_type[atom_type]

            # Assign scale factor and HCT radius using the atom element
            else:
                found = True

                # Assign scale factor
                if element.upper() in scale_by_element:
                    scale = scale_by_element[element.upper()]
                    scale = _check_bonds(scale, element.upper(),
                                         degree, parent)
                else:
                    found = False

                # Assign scale factor
                if element.upper() in radius_by_element:
                    radius = radius_by_element[element.upper()]
                else:
                    found = False

                # Returns scale and HCT radius if found, otherwise it raises a
                # warning and returns the default parameters
                if not found:
                    log = Logger()
                    log.warning('Warning: OBC parameters for '
                                + '{} {} '.format(atom_name, atom_type)
                                + 'NOT found in the template '
                                + 'database. Using default parameters')
                    radius, scale = float(0.80), float(2.0)

            return radius, scale

        from simtk import unit

        # Loop over atom types and names:
        radii = list()
        scales = list()
        for atom_name, atom_type in zip(OPLS_params['atom_names'],
                                        OPLS_params['atom_types']):
            atom_name = re.sub('_', ' ', atom_name)
            radius, scale = _find_GBSA_parameters_according_to(
                atom_name, atom_type,
                degree_by_name.get(atom_name),
                element_by_name.get(atom_name),
                parent_by_name.get(atom_name))

            radii.append(unit.Quantity(float(radius), unit.angstrom),)
            scales.append(scale)

        # Assign OBC parameters
        OPLS_params['GBSA_radii'] = radii
        OPLS_params['GBSA_scales'] = scales


class OpenFFOPLS2005ParameterWrapper(BaseParameterWrapper):
    """
    It defines a parameters wrapper for an hybrid OpenFF-OPLS2005 force
    field.
    """

    _name = 'Openff + OPLS2005'
