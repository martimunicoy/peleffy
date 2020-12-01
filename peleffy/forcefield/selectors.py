"""
This module handles the selection of the right force field or parameter
calculators, given a context.
"""


class ForceFieldSelector(object):
    """
    It defines a force field selector.
    """
    _FF_TYPES = {'OPLS2005': ('OPLS2005'),
                 'OpenFF': ('openff_unconstrained-1.3.0.offxml',
                            'openff_unconstrained-1.2.1.offxml',
                            'openff_unconstrained-1.2.0.offxml',
                            'openff_unconstrained-1.1.1.offxml',
                            'openff_unconstrained-1.1.0.offxml',
                            'openff_unconstrained-1.0.1.offxml',
                            'openff_unconstrained-1.0.0.offxml')}

    def get_by_name(self, forcefield_name):
        """
        Given a forcefield name, it returns the corresponding
        force field class.

        Parameters
        ----------
        forcefield_name : str
            The name of the requested forcefield

        Returns
        -------
        forcefield : a peleffy.forcefield.forcefield.ForceField
            The force field that corresponds to the supplied name

        Raises
        ------
        ValueError
            If the supplied forcefield_name is unknown

        Examples
        --------

        Use the force field selector to select a force field by name

        >>> from peleffy.forcefield import ForceFieldSelector

        >>> selector = ForceFieldSelector()
        >>> openff = selector.get_by_name('openff_unconstrained-1.2.1.offxml')

        """
        from peleffy.utils import Logger

        log = Logger()
        log.info(' - Loading \'{}\''.format(forcefield_name))

        from .forcefield import (OpenForceField, OPLS2005ForceField)

        if forcefield_name.upper() in self._FF_TYPES['OPLS2005']:
            return OPLS2005ForceField()
        elif forcefield_name in self._FF_TYPES['OpenFF']:
            return OpenForceField(forcefield_name=forcefield_name)
        else:
            raise ValueError('Invalid force field name')

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


class ChargeCalculatorSelector(object):
    """
    It defines a charge calculator selector.
    """

    from peleffy.forcefield.calculators import (OPLSChargeCalculator,
                                                Am1bccCalculator,
                                                GasteigerCalculator,
                                                DummyChargeCalculator)

    _AVAILABLE_TYPES = {'opls2005': OPLSChargeCalculator,
                        'am1bcc': Am1bccCalculator,
                        'gasteiger': GasteigerCalculator,
                        'dummy': DummyChargeCalculator
                        }

    def get_by_name(self, charge_method, molecule):
        """
        Given a charge method name, it returns the corresponding
        charge calculator class.

        Parameters
        ----------
        charge_method : str
            The name of the requested charge calculator
        molecule : a peleffy.topology.Molecule
            The peleffy's Molecule object whose partial charges will
            be calculated

        Returns
        -------
        charge_calculator : a PartialChargesCalculator
            The charge calculation method that will be employed to
            calculate partial charges

        Raises
        ------
        ValueError
            If the requested charge method is unknown
        """

        if charge_method.lower() not in self._AVAILABLE_TYPES:
            raise ValueError('Charge method \'{}\' '.format(charge_method)
                             + 'is unknown')

        charge_method = self._AVAILABLE_TYPES[charge_method.lower()]

        return charge_method(molecule)

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
