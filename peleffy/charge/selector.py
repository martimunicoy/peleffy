"""
This module handles the selection of the right charge calculator method.
"""


from peleffy.charge import (Am1bccCalculator, GasteigerCalculator,
                            OPLSChargeCalculator)


class ChargeCalculatorSelector(object):
    """
    It defines a charge calculator selector.
    """
    _AVAILABLE_TYPES = {'opls2005': OPLSChargeCalculator,
                        'am1bcc': Am1bccCalculator,
                        'gasteiger': GasteigerCalculator
                        }

    def get_by_name(self, charge_method):
        """
        Given a charge method name, it returns the corresponding
        charge calculator class.

        Parameters
        ----------
        charge_method : str
            The name of the requested charge calculator

        Returns
        -------
        charge_calculator : a PartialChargesCalculator object
            The charge calculation method that will be employed to calculate
            partial charges

        Raises
        ------
        ValueError
            If the requested charge method is unknown
        """

        if charge_method.lower() not in self._AVAILABLE_TYPES:
            raise ValueError('Charge method \'{}\' '.format(charge_method)
                             + 'is unknown')

        return self._AVAILABLE_TYPES[charge_method.lower()]

    def get_list(self):
        """
        It returns the complete list of available force fields.

        Returns
        -------
        forcefields : dict
            The complete list of available force fields grouped by
            force field type
        """
        return self._FF_TYPES
