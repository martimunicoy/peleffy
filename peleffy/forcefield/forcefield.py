"""
This module contains classes that handle force field representations.
"""


__all__ = ["OpenForceField", "OPLS2005ForceField",
           "OpenFFOPLS2005ForceField"]


class _BaseForceField(object):
    """
    It is the force field base class.
    """

    _type = ''

    def __init__(self, forcefield_name):
        """
        It initializes the base force field class.

        Parameters
        ----------
        forcefield_name : str
            The name of the forcefield
        """
        self._name = forcefield_name
        self._parameters = None

    def parameterize(self, molecule, force_parameterization=False):
        """
        It parameterizes the supplied molecule.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize
        force_parameterization : bool
            Whether to force a new parameterization instead of attempting
            to reuse parameters obtained in a previous parameterization,
            or not

        Returns
        -------
        parameters : an peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field

        """
        if self.parameters is None or force_parameterization:
            self._parameters = self._get_parameters(molecule)

        return self.parameters

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

    def _get_parameters(self, molecule):
        """
        It parameterizes the supplied molecule with the OpenFF toolkit.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
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

        from peleffy.forcefield import OpenForceFieldParameterWrapper

        return OpenForceFieldParameterWrapper.from_label_molecules(
            molecule, parameters)


class OPLS2005ForceField(_BaseForceField):
    """
    It defines an OPLS2005 force field.
    """

    _type = 'OPLS2005'

    def _get_parameters(self, molecule):
        """
        It parameterizes the supplied molecule with Schrodinger toolkit.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
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

        from peleffy.forcefield import OPLS2005ParameterWrapper

        return OPLS2005ParameterWrapper.from_ffld_output(molecule,
                                                         ffld_output)


class OpenFFOPLS2005ForceField(_BaseForceField):
    """
    It defines an hybrid force field from the combination of OpenFF and
    OPLS2005.
    """

    _type = 'OpenFF + OPLS2005'
    _selections = ('openff', 'opls2005')

    def __init__(self, forcefield_name):
        """
        It initializes the OpenFFOPLS2005ForceField class.

        Parameters
        ----------
        openff_name : str
            The name of the OpenFF forcefield to employ
        """
        self._openff = OpenForceField(forcefield_name)
        self._oplsff = OPLS2005ForceField('OPLS2005')
        super().__init__(self._openff.name + ' + ' + self._oplsff.name)

        # Set default selections
        self._nonbonding = 'openff'
        self._bonds = 'openff'
        self._angles = 'openff'
        self._torsions = 'openff'

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
        if selection.lower() not in self._selections:
            raise ValueError('Invalid selection: \'{}\'. '.format(selection)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        self._nonbonding = selection.lower()

    def set_bond_parameters(self, selection):
        """
        It sets which force field to use to assign parameters to bonds.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        if selection.lower() not in self._selections:
            raise ValueError('Invalid selection: \'{}\'. '.format(selection)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        self._bonds = selection.lower()

    def set_angle_parameters(self, selection):
        """
        It sets which force field to use to assign parameters to angles.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        if selection.lower() not in self._selections:
            raise ValueError('Invalid selection: \'{}\'. '.format(selection)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        self._angles = selection.lower()

    def set_torsion_parameters(self, selection):
        """
        It sets which force field to use to assign parameters to torsions.

        Parameters
        ----------
        selection : str
            The force field to employ. One of ('OpenFF', 'OPLS2005')
        """
        if selection.lower() not in self._selections:
            raise ValueError('Invalid selection: \'{}\'. '.format(selection)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        self._torsions = selection.lower()

    def _get_parameters(self, molecule):
        """
        It parameterizes the supplied molecule with both OpenFF and
        Schrodinger toolkits.

        Parameters
        ----------
        molecule : an peleffy.topology.Molecule
            The peleffy's Molecule object to parameterize

        Returns
        -------
        parameters : an OpenForceFieldParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field
        """
        from peleffy.forcefield import OpenFFOPLS2005ParameterWrapper

        hybrid_parameters = OpenFFOPLS2005ParameterWrapper()

        openff_parameters = self._openff.parameterize(molecule)
        oplsff_parameters = self._oplsff.parameterize(molecule)

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

        if self._torsions == 'openff':
            hybrid_parameters['propers'] = openff_parameters['propers']
            hybrid_parameters['impropers'] = openff_parameters['impropers']
        elif self._torsions == 'opls2005':
            hybrid_parameters['propers'] = oplsff_parameters['propers']
            hybrid_parameters['impropers'] = oplsff_parameters['impropers']
        else:
            raise ValueError('Invalid selection: '
                             + '\'{}\'. '.format(self._torsions)
                             + 'Valid selections are '
                             + '{}'.format(self._selections))

        return hybrid_parameters
