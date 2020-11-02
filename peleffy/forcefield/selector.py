"""
This module handles the selection of the right force field, given a context.
"""


class ForceFieldSelector(object):
    """
    It defines a force field selector.
    """
    _FF_TYPES = {'OPLS2005': ('OPLS2005'),
                 'OpenFF': ('openff_unconstrained-1.2.1.offxml',
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
        forcefield : an peleffy.forcefield.forcefield.ForceField object
            The force field object that corresponds to the supplied name

        Raises
        ------
        ValueError
            If the supplied forcefield_name is unknown
        """

        from .forcefield import (OpenForceField, OPLS2005ForceField)

        if forcefield_name.upper() in self._FF_TYPES['OPLS2005']:
            return OPLS2005ForceField(forcefield_name)
        elif forcefield_name in self._FF_TYPES['OpenFF']:
            return OpenForceField(forcefield_name)
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
