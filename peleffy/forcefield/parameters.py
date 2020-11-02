"""
This module provides classes to define force field parameters.
"""


__all__ = ["OpenForceFieldParameterWrapper",
           "OPLS2005ParameterWrapper",
           "OpenFFOPLS2005ParameterWrapper"]


from collections import defaultdict

from peleffy.utils import get_data_file_path


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

    def __init__(self, parameters_dict=dict()):
        """
        It initializes a parameter wrapper object.

        Parameters
        ----------
        parameters_dict : dict
             A dictionary keyed by parameter type that contains the
             parameters to store. Default is an empty dict
        """
        for key in self._keys:
            self[key] = list()

        for key, value in parameters_dict.items():
            self[key] = value

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


class OpenForceFieldParameterWrapper(BaseParameterWrapper):
    """
    It defines a parameters wrapper for an OpenFF force field.
    """

    _name = 'OpenFF'

    @staticmethod
    def from_label_molecules(molecule, parameters):
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
                     'spring_constant': k_by_idx[idxs] / 2.0,  # PELE works with half of the OFF's spring
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
                     'spring_constant': k_by_idx[idxs] / 2.0,  # PELE works with half of the OFF's spring
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

        return OpenForceFieldParameterWrapper(params)


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

                name_to_index[line[0:4]] = len(name_to_index)

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
                    {'atom1_idx': name_to_index[line[0:4]],
                     'atom2_idx': name_to_index[line[8:12]],
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
                    {'atom1_idx': name_to_index[line[0:4]],
                     'atom2_idx': name_to_index[line[8:12]],
                     'atom3_idx': name_to_index[line[16:20]],
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

                atom1_idx = name_to_index[line[0:4]]
                atom2_idx = name_to_index[line[8:12]]
                atom3_idx = name_to_index[line[16:20]]
                atom4_idx = name_to_index[line[24:28]]

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
                    {'atom1_idx': name_to_index[line[0:4]],
                     'atom2_idx': name_to_index[line[8:12]],
                     'atom3_idx': name_to_index[line[16:20]],
                     'atom4_idx': name_to_index[line[24:28]],
                     'periodicity': 2,
                     'phase': unit.Quantity(180.0, unit.degree),
                     'k': k / 2.0,  # PELE works with half of Schrodinger's force constant
                     'idivf': 1.0
                     })

        opls_parameters_wrapper = OPLS2005ParameterWrapper(params)
        OPLS2005ParameterWrapper._add_solvent_parameters(
            opls_parameters_wrapper)

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
    def _add_solvent_parameters(OPLS_params):
        """
        It add the solvent parameters to the OPLS parameters collection.

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


class OpenFFOPLS2005ParameterWrapper(BaseParameterWrapper):
    """
    It defines a parameters wrapper for an hybrid OpenFF-OPLS2005 force
    field.
    """

    _name = 'Openff + OPLS2005'
