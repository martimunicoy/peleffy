Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation

0.3.0 - Current development
-------------------------

This is still a preliminary version of the Open Force Field to PELE package. It includes some new features and improvements in documentation and unit testing.

New features
""""""""""""
- `PR #15 <https://github.com/martimunicoy/offpele/pull/15>`_: Adds a new method (Antechamber's gasteiger) to calculate partial charges.

Tests added
"""""""""""
- `PR #15 <https://github.com/martimunicoy/offpele/pull/15>`_: Adds tests ensuring that the run_offpele call from main and the partial charge calculators work as expected.


0.2.0
-----

This is a preliminary version of the Open Force Field to PELE package.

New features
""""""""""""

A first implementation of the package that allows to:

- Build a rotamer library for a small molecule using RDKit's API
- Build a template with the Molecular Mechanics' parameters for a small molecule using the Open Force Field Toolkit
- Assign the OBC implicit solvent parameters to a small molecule using the Open Force Field Toolkit
