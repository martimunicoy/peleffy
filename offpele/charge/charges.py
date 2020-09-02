"""
This module handles all classes and functions related with partial charge
calculators.
"""


from offpele.utils.toolkits import (AmberToolkitWrapper,
                                    ToolkitUnavailableException)


class _PartialChargesCalculator(object):
    """
    Base class for partial charges calculators.
    """

    _name = None
    _amber_toolkit = AmberToolkitWrapper()

    def __init__(self, molecule):
        """
        It initiates a PartialChargesCalculator object.

        Parameters
        ----------
        molecule : An offpele.topology.Molecule
            The partial charges of this Molecule object will be calculated
        """
        self._molecule = molecule

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


class Am1bccCalculator(_PartialChargesCalculator):
    """
    Implementation of the AM1-BCC partial charges calculator (using RDKit).
    """

    _name = 'am1bcc'


class GasteigerCalculator(_PartialChargesCalculator):
    """
    Implementation of the gasteiger partial charges calculator (using
    RDKit).
    """

    _name = 'gasteiger'


class OPLSChargeCalculator(_PartialChargesCalculator):
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
        charges : simtk.unit.Quantity
            The array of partial charges
        """

        try:
            OPLS_params = self.molecule.get_OPLS_parameters()
        except ToolkitUnavailableException:
            raise ToolkitUnavailableException(
                'OPLSChargeCalculator requires the Schrodinger '
                + 'Toolkit to obtain partial charges')

        return OPLS_params['charges']
