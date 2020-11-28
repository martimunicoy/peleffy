"""
This module contains classes that handle force field representations.
"""


__all__ = ["OpenForceField", "OPLS2005ForceField",
           "OpenFFOPLS2005ForceField"]


from peleffy.forcefield.selectors import ChargeCalculatorSelector


class _BaseForceField(object):
    """
    It is the force field base class.
    """

    _type = ''
    _default_charge_method = 'am1bcc'

    def __init__(self, forcefield_name):
        """
        It initializes the base force field class.

        Parameters
        ----------
        forcefield_name : str
            The name of the forcefield

        Examples
        --------

        Load a molecule and generate its parameters with an OpenFF
        force field

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')

        >>> from peleffy.forcefield import OpenForceField

        >>> openff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        >>> parameters = openff.parameterize(molecule)

        Load a molecule and generate its parameters with an OpenFF
        force field and 'gasteiger' partial charges

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')

        >>> from peleffy.forcefield import OpenForceField

        >>> openff = OpenForceField('openff_unconstrained-1.2.1.offxml')
        >>> parameters = openff.parameterize(molecule,
                                             charge_method='gasteiger')
        """
        self._name = forcefield_name

    def _get_charge_calculator(self, charge_method, molecule):
        """
        Given a charge method name, it returns the corresponding
        charge calculator class.

        Parameters
        ----------
        charge_method : str
            The name of the requested charge calculator
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize

        Returns
        -------
        charge_calculator : a PartialChargesCalculator object
            The charge calculation method that will be employed to calculate
            partial charges
        """

        # If charge method is not supplied, use force field's default
        if charge_method is None:
            charge_method = self._default_charge_method

        selector = ChargeCalculatorSelector()
        return selector.get_by_name(charge_method, molecule)

    def parameterize(self, molecule, charge_method=None):
        """
        It parameterizes the supplied molecule.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize
        charge_method : str
            The name of the charge method to employ

        Returns
        -------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field
        """

        # Assign parameters
        parameters = self._get_parameters(molecule)

        # Assign partial charges using the charge calculator object
        charge_calculator = self._get_charge_calculator(charge_method,
                                                        molecule)
        charge_calculator.assign_partial_charges(parameters)

        return parameters

    @property
    def type(self):
        """
        It returns the type of the force field.

        Returns
        -------
        type : str
            The type of the force field
        """
        return self._type

    @property
    def name(self):
        """
        It returns the name of the force field.

        Returns
        -------
        name : str
            The name of the force field
        """
        return self._name

    @property
    def parameters(self):
        """
        It returns the parameters obtained with the force field.

        Returns
        -------
        parameters : any class deriving from BaseParameterWrapper
            The parameters obtained with the force field
        """
        return self._parameters


class OpenForceField(_BaseForceField):
    """
    It defines an OpenFF force field.
    """

    _type = 'OpenFF'
    _default_charge_method = 'am1bcc'

    def _get_parameters(self, molecule):
        """
        It parameterizes the supplied molecule with the OpenFF toolkit.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize

        Returns
        -------
        parameters : an OpenForceFieldParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field
        """
        from peleffy.utils.toolkits import OpenForceFieldToolkitWrapper

        openforcefield_toolkit = OpenForceFieldToolkitWrapper()
        parameters = openforcefield_toolkit.get_parameters_from_forcefield(
            self.name, molecule)

        from peleffy.forcefield.parameters \
            import OpenForceFieldParameterWrapper

        return OpenForceFieldParameterWrapper.from_label_molecules(
            molecule, parameters, self.name)


class OPLS2005ForceField(_BaseForceField):
    """
    It defines an OPLS2005 force field.
    """

    _type = 'OPLS2005'
    _default_charge_method = 'opls2005'

    def __init__(self):
        """
        It initializes the OPLS2005 force field class.

        Examples
        --------

        Load a molecule and generate its parameters with the OPLS2005
        force field

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')

        >>> from peleffy.forcefield import OPLS2005ForceField

        >>> opls2005 = OPLS2005ForceField()
        >>> parameters = opls2005.parameterize(molecule)

        """

        super().__init__(self._type)

    def _get_parameters(self, molecule):
        """
        It parameterizes the supplied molecule with Schrodinger toolkit.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize

        Returns
        -------
        parameters : an OpenForceFieldParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field
        """
        from peleffy.utils.toolkits import SchrodingerToolkitWrapper

        schrodinger_toolkit = SchrodingerToolkitWrapper()
        ffld_output = schrodinger_toolkit.run_ffld_server(molecule)

        from peleffy.forcefield.parameters \
            import OPLS2005ParameterWrapper

        return OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                         ffld_output)


class OpenFFOPLS2005ForceField(_BaseForceField):
    """
    It defines an hybrid force field from the combination of OpenFF and
    OPLS2005.
    """

    _type = 'OpenFF + OPLS2005'
    _selections = ('openff', 'opls2005')
    _default_charge_method = 'am1bcc'

    def __init__(self, forcefield_name):
        """
        It initializes the OpenFFOPLS2005ForceField class.

        Parameters
        ----------
        openff_name : str
            The name of the OpenFF forcefield to employ

        Examples
        --------

        Load a molecule and generate its parameters with an hybrid
        OpenFF-OPLS2005 force field

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')

        >>> from peleffy.forcefield import OpenFFOPLS2005ForceField

        >>> hybridff = OpenFFOPLS2005ForceField('openff_unconstrained-1.2.1.offxml')
        >>> parameters = hybridff.parameterize(molecule)

        We can customized it by selecting which force field employ for each
        parameter type. For example, below we are parameterizing all
        force field elements with OPLS2005 except for the dihedrals

        >>> from peleffy.topology import Molecule

        >>> molecule = Molecule('molecule.pdb')

        >>> from peleffy.forcefield import OpenFFOPLS2005ForceField

        >>> hybridff = OpenFFOPLS2005ForceField('openff_unconstrained-1.2.1.offxml')

        >>> hybridff.set_nonbonding_parameters('OPLS2005')
        >>> hybridff.set_bond_parameters('OPLS2005')
        >>> hybridff.set_angle_parameters('OPLS2005')
        >>> hybridff.set_dihedral_parameters('OpenFF')

        >>> parameters = hybridff.parameterize(molecule)

        """
        self._openff = OpenForceField(forcefield_name)
        self._oplsff = OPLS2005ForceField()
        super().__init__(self._openff.name + ' + ' + self._oplsff.name)

        # Set default selections
        self._nonbonding = 'openff'
        self._bonds = 'openff'
        self._angles = 'openff'
        self._dihedrals = 'openff'

    def _check_selection(self, selection):
        """
        It performs a safety check on the parameter selection.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')

        Raises
        ------
        ValueError if the selection string is unknown
        """
        if selection.lower() not in self._selections:
            raise ValueError('Invalid selection: \'{}\'. '.format(selection)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

    def set_nonbonding_parameters(self, selection):
        """
        It sets which force field to use to assign nonbonding parameters.
        Nonbonding parameters are: 'atom_types', 'sigmas', 'epsilons',
        'SGB_radii', 'vdW_radii', 'gammas', 'alphas', 'GBSA_radii' and
        'GBSA_scales'.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        self._check_selection(selection)

        self._nonbonding = selection.lower()

    def set_bond_parameters(self, selection):
        """
        It sets which force field to use to assign parameters to bonds.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        self._check_selection(selection)

        self._bonds = selection.lower()

    def set_angle_parameters(self, selection):
        """
        It sets which force field to use to assign parameters to angles.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        self._check_selection(selection)

        self._angles = selection.lower()

    def set_dihedral_parameters(self, selection):
        """
        It sets which force field to use to assign parameters to dihedrals.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        self._check_selection(selection)

        self._dihedrals = selection.lower()

    def _get_parameters(self, molecule):
        """
        It parameterizes the supplied molecule with both OpenFF and
        Schrodinger toolkits.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize

        Returns
        -------
        parameters : an OpenForceFieldParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field
        """
        from peleffy.forcefield.parameters \
            import OpenFFOPLS2005ParameterWrapper

        hybrid_parameters = OpenFFOPLS2005ParameterWrapper(
            forcefield_name='{} + {}'.format(self._openff.name,
                                             self._oplsff.name))

        openff_parameters = self._openff.parameterize(molecule,
                                                      charge_method='dummy')
        print(openff_parameters['sigmas'])
        oplsff_parameters = self._oplsff.parameterize(molecule)
        print(oplsff_parameters['sigmas'])

        if self._nonbonding == 'openff':
            hybrid_parameters['atom_names'] = openff_parameters['atom_names']
            hybrid_parameters['atom_types'] = openff_parameters['atom_types']
            hybrid_parameters['charges'] = openff_parameters['charges']
            hybrid_parameters['sigmas'] = openff_parameters['sigmas']
            hybrid_parameters['epsilons'] = openff_parameters['epsilons']
            hybrid_parameters['SGB_radii'] = openff_parameters['SGB_radii']
            hybrid_parameters['vdW_radii'] = openff_parameters['vdW_radii']
            hybrid_parameters['gammas'] = openff_parameters['gammas']
            hybrid_parameters['alphas'] = openff_parameters['alphas']
            hybrid_parameters['GBSA_radii'] = openff_parameters['GBSA_radii']
            hybrid_parameters['GBSA_scales'] = openff_parameters['GBSA_scales']
        elif self._nonbonding == 'opls2005':
            hybrid_parameters['atom_names'] = oplsff_parameters['atom_names']
            hybrid_parameters['atom_types'] = oplsff_parameters['atom_types']
            hybrid_parameters['charges'] = oplsff_parameters['charges']
            hybrid_parameters['sigmas'] = oplsff_parameters['sigmas']
            hybrid_parameters['epsilons'] = oplsff_parameters['epsilons']
            hybrid_parameters['SGB_radii'] = oplsff_parameters['SGB_radii']
            hybrid_parameters['vdW_radii'] = oplsff_parameters['vdW_radii']
            hybrid_parameters['gammas'] = oplsff_parameters['gammas']
            hybrid_parameters['alphas'] = oplsff_parameters['alphas']
            hybrid_parameters['GBSA_radii'] = oplsff_parameters['GBSA_radii']
            hybrid_parameters['GBSA_scales'] = oplsff_parameters['GBSA_scales']
        else:
            raise ValueError('Invalid selection: '
                             + '\'{}\'. '.format(self._nonbonding)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        if self._bonds == 'openff':
            hybrid_parameters['bonds'] = openff_parameters['bonds']
        elif self._bonds == 'opls2005':
            hybrid_parameters['bonds'] = oplsff_parameters['bonds']
        else:
            raise ValueError('Invalid selection: '
                             + '\'{}\'. '.format(self._bonds)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        if self._angles == 'openff':
            hybrid_parameters['angles'] = openff_parameters['angles']
        elif self._angles == 'opls2005':
            hybrid_parameters['angles'] = oplsff_parameters['angles']
        else:
            raise ValueError('Invalid selection: '
                             + '\'{}\'. '.format(self._angles)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        if self._dihedrals == 'openff':
            hybrid_parameters['propers'] = openff_parameters['propers']
            hybrid_parameters['impropers'] = openff_parameters['impropers']
        elif self._dihedrals == 'opls2005':
            hybrid_parameters['propers'] = oplsff_parameters['propers']
            hybrid_parameters['impropers'] = oplsff_parameters['impropers']
        else:
            raise ValueError('Invalid selection: '
                             + '\'{}\'. '.format(self._dihedrals)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        return hybrid_parameters
