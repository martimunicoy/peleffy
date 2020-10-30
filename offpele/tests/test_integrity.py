"""
This module contains integrity tests that check that all offpele's modules
are available.
"""

import pytest


class TestIntegrity(object):
    """Integrity test."""

    def test_modules(self):
        try:
            import offpele
            from offpele import charge
            from offpele import forcefield
            from offpele import solvent
            from offpele import template
            from offpele import topology
            from offpele import utils

        except ImportError as e:
            raise AssertionError('The following offpele module is missing: '
                                 + str(e))
