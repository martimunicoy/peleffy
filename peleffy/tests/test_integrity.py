"""
This module contains integrity tests that check that all peleffy's modules
are available.
"""

import pytest


class TestIntegrity(object):
    """Integrity test."""

    def test_modules(self):
        try:
            import peleffy
            from peleffy import charge
            from peleffy import forcefield
            from peleffy import solvent
            from peleffy import template
            from peleffy import topology
            from peleffy import utils

        except ImportError as e:
            raise AssertionError('The following peleffy module is missing: '
                                 + str(e))
