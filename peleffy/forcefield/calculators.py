"""
This module contains different parameter calculators for force fields.
"""


from simtk import unit

from peleffy.utils.toolkits import ToolkitUnavailableException


class _PartialChargeCalculator(object):
    """
    Base class for partial charge calculators.
    """

    from peleffy.utils.toolkits import AmberToolkitWrapper

    _name = None
    _amber_toolkit = AmberToolkitWrapper()

    def __init__(self, molecule):
        """
        It initiates a PartialChargeCalculator object.

        Parameters
        ----------
        molecule : a peleffy.topology.Molecule
            The partial charges of this Molecule object will be calculated
        """
        self._molecule = molecule

    @property
    def molecule(self):
        """
        The peleffy's Molecule.

        Returns
        -------
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object
        """
        return self._molecule

    @property
    def name(self):
        return self._name

    def get_partial_charges(self):
        """
        It returns the partial charges that correspond to the molecule's
        atoms.

        Returns
        -------
        charges : simtk.unit.Quantity
            The array of partial charges
        """

        return self._amber_toolkit.compute_partial_charges(self.molecule,
                                                           method=self.name)

    def assign_partial_charges(self, parameters):
        """
        Given a parameter wrapper, it computes the partial charges and
        assigns them to it.

        Parameters
        ----------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field and where the partial charges
            will be assigned
        """

        partial_charges = self.get_partial_charges()
        parameters['charges'] = partial_charges


class Am1bccCalculator(_PartialChargeCalculator):
    """
    Implementation of the AM1-BCC partial charge calculator (using RDKit).
    """

    _name = 'am1bcc'


class GasteigerCalculator(_PartialChargeCalculator):
    """
    Implementation of the gasteiger partial charge calculator (using
    RDKit).
    """

    _name = 'gasteiger'


class OPLSChargeCalculator(_PartialChargeCalculator):
    """
    Implementation of the calculator of OPLS partial charges (using
    Schrodinger's ffld_server)
    """

    _name = 'OPLS2005'

    def get_partial_charges(self):
        """
        It returns the partial charges that correspond to the molecule's
        atoms.

        Returns
        -------
        partial_charges : simtk.unit.Quantity
            The array of partial charges
        """

        from peleffy.utils.toolkits import SchrodingerToolkitWrapper

        try:
            schrodinger_toolkit_wrapper = SchrodingerToolkitWrapper()
            ffld_output = schrodinger_toolkit_wrapper.run_ffld_server(
                self.molecule)
        except ToolkitUnavailableException:
            raise ToolkitUnavailableException(
                'OPLSChargeCalculator requires the Schrodinger '
                + 'Toolkit to obtain partial charges')

        from peleffy.forcefield.parameters import OPLS2005ParameterWrapper
        parameters = OPLS2005ParameterWrapper.from_ffld_output(
            self.molecule, ffld_output)

        return parameters['charges']

    def assign_partial_charges(self, parameters):
        """
        Given a parameter wrapper, it computes the partial charges and
        assigns them to it.

        Parameters
        ----------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field and where the partial charges
            will be assigned
        """
        import peleffy

        if not isinstance(
                parameters,
                peleffy.forcefield.parameters.OPLS2005ParameterWrapper):
            partial_charges = self.get_partial_charges()
            parameters['charges'] = partial_charges


class DummyChargeCalculator(_PartialChargeCalculator):
    """
    Implementation of a dummy charge calculator that will not perform
    any partial charge calculation.
    """

    _name = 'OPLS2005'

    def get_partial_charges(self):
        """
        It returns the partial charges that correspond to the molecule's
        atoms.

        Returns
        -------
        partial_charges : simtk.unit.Quantity
            The array of partial charges
        """

        # Extract number of atoms in molecule
        n_atoms = len(self.molecule.get_pdb_atom_names())

        return [unit.Quantity(0.0, unit.elementary_charge), ] * n_atoms

    def assign_partial_charges(self, parameters):
        """
        Given a parameter wrapper, it computes the partial charges and
        assigns them to it.

        Parameters
        ----------
        parameters : a peleffy.forcefield.parameters.BaseParameterWrapper object
            The parameter wrapper containing the parameters generated
            with the current force field and where the partial charges
            will be assigned
        """
        import peleffy

        if not isinstance(
                parameters,
                peleffy.forcefield.parameters.OPLS2005ParameterWrapper):
            partial_charges = self.get_partial_charges()
            parameters['charges'] = partial_charges
