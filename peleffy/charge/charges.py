"""
This module handles all classes and functions related with partial charge
calculators.
"""


from peleffy.utils.toolkits import (AmberToolkitWrapper,
                                    ToolkitUnavailableException)


class _PartialChargeCalculator(object):
    """
    Base class for partial charge calculators.
    """

    _name = None
    _amber_toolkit = AmberToolkitWrapper()

    def __init__(self, molecule):
        """
        It initiates a PartialChargeCalculator object.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            The partial charges of this Molecule object will be calculated
        """
        self._molecule = molecule

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

    _name = 'OPLS'

    def get_partial_charges(self):
        """
        It returns the partial charges that correspond to the molecule's
        atoms.

        Returns
        -------
        partial_charges : simtk.unit.Quantity
            The array of partial charges
        """

        parameters = self.molecule.parameters

        # Only attempt to compute partial charges if the pure OpenFF has
        # been employed. Otherwise, OPLS2005 partial charges will already
        # be calculated and stored to molecule's parameters
        if parameters is None or parameters.name == 'OpenFF':
            try:
                from peleffy.utils.toolkits import SchrodingerToolkitWrapper

                schrodinger_toolkit_wrapper = SchrodingerToolkitWrapper()
                ffld_output = schrodinger_toolkit_wrapper.run_ffld_server(
                    self.molecule)
            except ToolkitUnavailableException:
                raise ToolkitUnavailableException(
                    'OPLSChargeCalculator requires the Schrodinger '
                    + 'Toolkit to obtain partial charges')

            from peleffy.forcefield import OPLS2005ParameterWrapper
            parameters = OPLS2005ParameterWrapper.from_ffld_output(
                self.molecule, ffld_output)

        """
        partial_charges = list()
        for partial_charge in parameters['charges']:
            value = partial_charge.value_in_unit(unit.elementary_charge)
            partial_charges.append(value)


        partial_charges = unit.Quantity(np.array(partial_charges),
                                        unit.elementary_charge)
        """

        return parameters['charges']
