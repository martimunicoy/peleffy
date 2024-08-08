Release History
===============

Releases follow the ``major.minor.micro`` scheme recommended by `PEP440 <https://www.python.org/dev/peps/pep-0440/#final-releases>`_, where

* ``major`` increments denote a change that may break API compatibility with previous ``major`` releases
* ``minor`` increments add features but do not break API compatibility
* ``micro`` increments represent bugfix releases or improvements in documentation


1.5.1 - Important bug fix for alchemistry
-----------------------------------------

This is a minor release of peleffy that corrects an important problem in alchemical
templates.

Bugfixes
""""""""
- `PR #183 <https://github.com/martimunicoy/peleffy/pull/183>`_: bug fixes for alchemical topologies


1.5.0 - OpenFF upgrade and enhancements for alchemistry
-------------------------------------------------------

This is a minor release of peleffy that supports the latest version of openff-toolkit 0.10. It also
contains additional enhancements for the alchemistry module.

New features
""""""""""""
- `PR #182 <https://github.com/martimunicoy/peleffy/pull/182>`_: support openff-toolkit 0.10.7 and openff-forcefield 2.1.0
- `PR #182 <https://github.com/martimunicoy/peleffy/pull/182>`_: better alchemical structure alignment
- `PR #182 <https://github.com/martimunicoy/peleffy/pull/182>`_: more lambda types available (vdw1 and vdw2)

Tests added
"""""""""""
- `PR #182 <https://github.com/martimunicoy/peleffy/pull/182>`_: corrections for alchemical and mapper tests
- `PR #182 <https://github.com/martimunicoy/peleffy/pull/182>`_: new test to check forcefield names


1.4.5 - Corrections for alchemistry
-----------------------------------

This is a micro release of peleffy that contains corrections for alchemistry module.

Bugfixes
""""""""
- `PR #180 <https://github.com/martimunicoy/peleffy/pull/180>`_: bug fixes for alchemical topologies


1.4.4 - Corrections for alchemistry and new charge calculator
-------------------------------------------------------------

This is a micro release of peleffy that contains corrections for alchemistry module. A new charge calculator
called Mulliken is also available.

New features
""""""""""""
- `PR #174 <https://github.com/martimunicoy/peleffy/pull/174>`_: minor fix in the logger to prevent conflicts with external loggers
- `PR #176 <https://github.com/martimunicoy/peleffy/pull/176>`_: new charge calculator called Mulliken
- `PR #177 <https://github.com/martimunicoy/peleffy/pull/177>`_: new method to save alchemical mapping as a PNG file

Bugfixes
""""""""
- `PR #177 <https://github.com/martimunicoy/peleffy/pull/177>`_: bug fixes for alchemical solvent templates and affected tests

Tests added
"""""""""""
- `PR #176 <https://github.com/martimunicoy/peleffy/pull/176>`_: new tests to validate the new charge calculator


1.4.3 - Minor improvements for CLI arguments and ffld_server
------------------------------------------------------------

This is a micro release of peleffy that contains minor improvements for the CLI interface and the usage of the ffld_server.

New features
""""""""""""
- `PR #171 <https://github.com/martimunicoy/peleffy/pull/171>`_: improvements for the CLI interface. Also, ffld_server will not raise an exception but any warning or error found will be raised by peleffy's logger


1.4.2 - Compatibility fixes for latest RDKit and Schrodinger versions
---------------------------------------------------------------------

This is a micro release of peleffy that fixes several bugs related with latest versions of RDKit and Schrodinger. Affected modules are rotamer libraries, the runner and parser of the ffld server, and the alchemistry package. It also adds some minor improvements to log handlers.

New features
""""""""""""
- `PR #166 <https://github.com/martimunicoy/peleffy/pull/166>`_: New options for log handlers.

Bugfixes
""""""""
- `PR #167 <https://github.com/martimunicoy/peleffy/pull/167>`_: Bug fixes for rotamer libraries and affected tests
- `PR #168 <https://github.com/martimunicoy/peleffy/pull/168>`_: Compatibility changes for the ffld server shipped with latest Schrodinger version
- `PR #169 <https://github.com/martimunicoy/peleffy/pull/169>`_: Support new RDKit versions.


1.4.1 - Bug fixes for heteromolecules extraction
---------------------------------------------------------

This is a micro release of peleffy that fixes a bug that affected heteromolecules extraction.

Bugfixes
""""""""
- `PR #163 <https://github.com/martimunicoy/peleffy/pull/163>`_: Fix problem when extracting heteromolecules from a specific chain

Tests added
"""""""""""
- `PR #163 <https://github.com/martimunicoy/peleffy/pull/163>`_: Adds a test to validate heteromolecule extraction


1.4.0 - Alchemistry and AMBER support
---------------------------------------------------------

This is a minor release of peleffy that adds a new module to generate alchemical templates. It also adds support for AMBER's implementation in PELE.

New features
""""""""""""
- `PR #154 <https://github.com/martimunicoy/peleffy/pull/154>`_: Introduces Alchemizer module to generate hybrid topologies
- `PR #155 <https://github.com/martimunicoy/peleffy/pull/155>`_: Adds support for PELE's AMBER with a new Impact template
- `PR #162 <https://github.com/martimunicoy/peleffy/pull/162>`_: Upgrade to openff-toolkit 0.10.1

Bugfixes
""""""""
- `PR #158 <https://github.com/martimunicoy/peleffy/pull/158>`_: Fix minor bug when using the --chain flag and introduces checks for the input PDB in the peleffy.main module
- `PR #159 <https://github.com/martimunicoy/peleffy/pull/159>`_: Fix issues with long atom numbers and heteromolecules extraction

Tests added
"""""""""""
- `PR #154 <https://github.com/martimunicoy/peleffy/pull/154>`_: Adds a collection of tests for Alchemizer module
- `PR #155 <https://github.com/martimunicoy/peleffy/pull/155>`_: Extends the tests for utils module and introduces new tests for the new AMBER-compatible Impact template
- `PR #158 <https://github.com/martimunicoy/peleffy/pull/158>`_: Extends the tests for the new checks in the peleffy.main module


1.3.4 - OpenFF-2.0 Support
---------------------------------------------------------

This is a micro release of peleffy that adds support for the new openff-2.0.0. It also solves minor bugs in the OPLS2005 parametrization.

New features
""""""""""""
- `PR #151 <https://github.com/martimunicoy/peleffy/pull/151>`_: Add support for openff-2.0.0.

Bugfixes
""""""""
- `PR #149 <https://github.com/martimunicoy/peleffy/pull/149>`_: Minor error when parsing ffld output file.
- `PR #153 <https://github.com/martimunicoy/peleffy/pull/153>`_: Fix parameters inconsistencies.


1.3.3 - Explicit hydrogens support
---------------------------------------------------------

This is a micro release of peleffy that includes support for the new OpenFF flag to manage explicit and implicit hydrogen atoms.

New features
""""""""""""
- `PR #146 <https://github.com/martimunicoy/peleffy/pull/146>`_: Adds support for the new explicit hydrogens flag

Tests added
"""""""""""
- `PR #146 <https://github.com/martimunicoy/peleffy/pull/146>`_: New test to check new explicit hydrogens flag


1.3.2 - Migration and support for openff.toolkit
---------------------------------------------------------

This is a micro release of peleffy that includes a migration to openff.toolkit to support future releases.

New features
""""""""""""
- `PR #144 <https://github.com/martimunicoy/peleffy/pull/144>`_: Migrated openforcefield imports to openff.toolkit


1.3.1 - PELE Platform support
-----------------------------

This is a micro release of peleffy that includes minor adjustments for the PELE Platform and other small fixes.

New features
""""""""""""
- `PR #142 <https://github.com/martimunicoy/peleffy/pull/142>`_: Minor adjustments to facilitate platform compatibility.

Bugfixes
""""""""
- `PR #143 <https://github.com/martimunicoy/peleffy/pull/143>`_: Minor error when parsing Impact templates.


1.3.0 - BCE conformations and automatic heteromolecules extraction
------------------------------------------------------------------

This is a minor release of peleffy that includes a new method to read ligand conformations from the BCE server and prepare the input files for PELE. It also contains a new PDBFile class that allows the user to automatically load all the heteromolecules from a PDB file.

New features
""""""""""""
- `PR #135 <https://github.com/martimunicoy/peleffy/pull/135>`_: New class to load conformations from the BCE server (Bioactive Conformational Ensemble) and generate the required input files for PELE.
- `PR #137 <https://github.com/martimunicoy/peleffy/pull/137>`_: New PDB class that allows to handle an input PDB file with multiple heteromolecules. 

Tests added
"""""""""""
- `PR #135 <https://github.com/martimunicoy/peleffy/pull/135>`_: Adds tests to validate the new BCEConformations class.
- `PR #137 <https://github.com/martimunicoy/peleffy/pull/137>`_: Adds tests to validate the new PDBFile class.


1.2.1 - API Documentation and improvements
------------------------------------------

This is a micro release of peleffy that includes and new method to load parameters from a JSON file and solves different bugs in the documentation and the OPLS parametrization. 

New features
""""""""""""
- `PR #131 <https://github.com/martimunicoy/peleffy/pull/131>`_: New method to load parameters from a JSON file.

Tests added
"""""""""""
- `PR #131 <https://github.com/martimunicoy/peleffy/pull/131>`_: Adds tests to validate the new method to load parameters from a JSON file.

Bugfixes
""""""""
- `PR #129 <https://github.com/martimunicoy/peleffy/pull/129>`_: Some format errors in the API documentation are fixed. Links to the PELE documentation are updated.
- `PR #134 <https://github.com/martimunicoy/peleffy/pull/134>`_: Fixes bug when parsing the parameters of the ligand when OPLS is used to parameterize. 


1.2.0 - New tools for parameters and templates
----------------------------------------------

This is a minor release of peleffy that includes new useful tools to handle parameters and their templates more easily. It also supports the newest version of the OpenForceField toolkit, which is 0.8.3.

New features
""""""""""""
- `PR #117 <https://github.com/martimunicoy/peleffy/pull/117>`_: New method to assign external partial charges.
- `PR #118 <https://github.com/martimunicoy/peleffy/pull/118>`_: New method to load parameters from an Impact Template.
- `PR #119 <https://github.com/martimunicoy/peleffy/pull/119>`_: Adds explanatory error message when using an invalid Impact Template in the from_impact_template method.
- `PR #119 <https://github.com/martimunicoy/peleffy/pull/119>`_: Supports Openforcefield-0.8.3 .
- `PR #126 <https://github.com/martimunicoy/peleffy/pull/126>`_: Allows the Solvent class to be compatible with multiple topologies. 

Bugfixes
""""""""
- `PR #125 <https://github.com/martimunicoy/peleffy/pull/125>`_: A bad index slicing in the molecule.Molecule._pdb_checkup() is now fixed. 

Tests added
"""""""""""
- `PR #117 <https://github.com/martimunicoy/peleffy/pull/117>`_: Adds tests to validate the MAE parse for external partial charges.
- `PR #118 <https://github.com/martimunicoy/peleffy/pull/118>`_: Adds tests to validate the new method to load parameters from an Impact Template.
- `PR #119 <https://github.com/martimunicoy/peleffy/pull/119>`_: Adds tests for the new error message when using an invalid Impact Template in the from_impact_template method.
- `PR #126 <https://github.com/martimunicoy/peleffy/pull/126>`_: Adds tests for the new compatibility of the Solvent class with multiple topologies. 


1.1.0 - Improvements in parameterization API, OBC template for OPLS2005 and Molecule initializators
---------------------------------------------------------------------------------------------------

This minor release introduces improvements to the parameterization API of peleffy. It also integrates the parameterization of OBC radii and scale factors required by the OPLS2005 implementation of PELE. It also improves the initialization of the Molecule class with a PDB checking and fixer and taking RDKit and OpenFF molecular representations as input. It also adds support for the new openff-1.3.0. Besides, it fixes a serious bug in the atom ordering of the Impact template that affected PELE's side chain prediction algorithm.

New features
""""""""""""
- `PR #86 <https://github.com/martimunicoy/peleffy/pull/86>`_: New method to check the input PDB prior building the molecule.
- `PR #88 <https://github.com/martimunicoy/peleffy/pull/88>`_: New method to retrieve atom degrees with RDKit.
- `PR #90 <https://github.com/martimunicoy/peleffy/pull/90>`_: Add support for openff-1.3.0.
- `PR #92 <https://github.com/martimunicoy/peleffy/pull/92>`_: New parameter to skip the stereochemistry assignment (and the checking from the OpenFF toolkit).
- `PR #94 <https://github.com/martimunicoy/peleffy/pull/94>`_: New method for the OPLS OBC parameters.
- `PR #100 <https://github.com/martimunicoy/peleffy/pull/100>`_: New writer for the OPLS OBC parameters.
- `PR #106 <https://github.com/martimunicoy/peleffy/pull/106>`_: New method to initialize a Molecule object directly from an RDKit and OpenFF molecular representations.
- `PR #112 <https://github.com/martimunicoy/peleffy/pull/112>`_: New method to fix an input PDB file with no atomic element identifiers.

Bugfixes
""""""""
- `PR #107 <https://github.com/martimunicoy/peleffy/pull/107>`_: A bad ordering of the atoms in the Impact template generated by peleffy is now fixed.
- `PR #115 <https://github.com/martimunicoy/peleffy/pull/115>`_: The rotamer library writer now writes molecule's tag not its name, as expected.

API-breaking changes
""""""""""""""""""""
- `PR #94 <https://github.com/martimunicoy/peleffy/pull/94>`_: Methods to write to a file are given a unique standard name, to_file(), to simplify the API.
- `PR #97 <https://github.com/martimunicoy/peleffy/pull/97>`_: The parameterization API changes and a new Topology class is used as a container for all the topological elements.
- `PR #103 <https://github.com/martimunicoy/peleffy/pull/103>`_: The OpenFF output of PELE changes from DataLocal/Templates/OFF/Parsley/HeteroAtoms/ to DataLocal/Templates/OpenFF/Parsley/.

Tests added
"""""""""""
- `PR #88 <https://github.com/martimunicoy/peleffy/pull/88>`_: Adds tests to validate the atom degrees getter.
- `PR #86 <https://github.com/martimunicoy/peleffy/pull/86>`_: Adds tests to validate the PDB check up.
- `PR #90 <https://github.com/martimunicoy/peleffy/pull/90>`_: General validation of supported force fields.
- `PR #92 <https://github.com/martimunicoy/peleffy/pull/92>`_: New test to check the behaviour of the allow_undefined_stereo parameter.
- `PR #94 <https://github.com/martimunicoy/peleffy/pull/94>`_: Adds tests to validate the OPLS OBC parameters generator.
- `PR #97 <https://github.com/martimunicoy/peleffy/pull/97>`_: Includes tests for the new Topology container class.
- `PR #100 <https://github.com/martimunicoy/peleffy/pull/100>`_: Adds tests to validate the solvent template writers.
- `PR #106 <https://github.com/martimunicoy/peleffy/pull/106>`_: Adds tests to check the RDKit and OpenFF molecular initializers.
- `PR #112 <https://github.com/martimunicoy/peleffy/pull/112>`_: Adds one test to check the new PDB fixer method.


1.0.0 - Full compatibility for OPLS2005 and OpenFF dihedrals and rotamer library improvements
---------------------------------------------------------------------------------------------

This major release renames the package to peleffy as it now supports both OpenFF and OPLS2005 force fields. Therefore, this release extends the compatibility of peleffy to fully handle OPLS2005 parameters. Some unsupported OpenFF dihedrals now can be handled. Besides, it includes a functionality to generate rotamer libraries with core constraints to allow the user to force an atom to be in the core.

New features
""""""""""""
- `PR #56 <https://github.com/martimunicoy/peleffy/pull/56>`_: Dynamic output path handler.
- `PR #62 <https://github.com/martimunicoy/peleffy/pull/62>`_: New functionality to generate rotamer libraries forcing an atom to be in the core.
- `PR #63 <https://github.com/martimunicoy/peleffy/pull/63>`_: Enhancements to the core constraints to allow the selection of multiple core atoms.
- `PR #66 <https://github.com/martimunicoy/peleffy/pull/66>`_: Full compatibility with OpenFF dihedrals.
- `PR #69 <https://github.com/martimunicoy/peleffy/pull/69>`_: Full compatibility with OPLS2005 force field.
- `PR #85 <https://github.com/martimunicoy/peleffy/pull/85>`_: Package is renamed to peleffy.

Bugfixes
""""""""
- `PR #74 <https://github.com/martimunicoy/peleffy/pull/74>`_: Corrects wrong assignment of PDB atom names when using the OPLS2005 force field.
- `PR #79 <https://github.com/martimunicoy/peleffy/pull/79>`_: Corrects error with missing modules in the Conda installation.
- `PR #82 <https://github.com/martimunicoy/peleffy/pull/82>`_: Corrects a bug that caused some important propers obtained with OPLS2005 to be missing.
- `PR #84 <https://github.com/martimunicoy/peleffy/pull/84>`_: Corrects a bug that caused unparameterized Molecules to be undetected.

Tests added
"""""""""""
- `PR #56 <https://github.com/martimunicoy/peleffy/pull/56>`_: Adds tests to validate the new output path handler.
- `PR #62 <https://github.com/martimunicoy/peleffy/pull/62>`_: Adds tests to validate the new rotamer library with core constraints.
- `PR #63 <https://github.com/martimunicoy/peleffy/pull/63>`_: More tests are added to validate the new rotamer library with core constraints.
- `PR #66 <https://github.com/martimunicoy/peleffy/pull/66>`_: Adds tests to validate the handling of non standard dihedrals.
- `PR #69 <https://github.com/martimunicoy/peleffy/pull/69>`_: Adds tests to validate the integration of OPLS2005 force field.
- `PR #70 <https://github.com/martimunicoy/peleffy/pull/70>`_: Adds tests to validate main CLI.
- `PR #84 <https://github.com/martimunicoy/peleffy/pull/840>`_: Adds tests to validate the Impact class.


0.3.1 - General stability improvements
--------------------------------------

This is a micro release that includes general bug fixes and stability improvements. It is still a preliminary version of the Open Force Field to PELE package which is under development.

New features
""""""""""""
- `PR #52 <https://github.com/martimunicoy/peleffy/pull/52>`_: Molecule connectivity can be assigned from an RDKit molecular template when loading it from a PDB file without connectivity.
- `PR #55 <https://github.com/martimunicoy/peleffy/pull/55>`_: Standard output prints follow the logging hierarchy and can be modified by the user.
- `PR #59 <https://github.com/martimunicoy/peleffy/pull/59>`_: Set alternative conformers to the peleffy's molecule representation.

Bugfixes
""""""""
- `PR #48 <https://github.com/martimunicoy/peleffy/pull/48>`_: Fixes CLI's default output paths.
- `PR #58 <https://github.com/martimunicoy/peleffy/pull/58>`_: Fixes unconsistencies between PDB residue name and molecule tag.

Tests added
"""""""""""
- `PR #48 <https://github.com/martimunicoy/peleffy/pull/48>`_: Adds tests to validate the assignment of the default output paths.
- `PR #52 <https://github.com/martimunicoy/peleffy/pull/52>`_: Adds tests to validate the initialization using a connectivity template.
- `PR #55 <https://github.com/martimunicoy/peleffy/pull/55>`_: Adds tests for the new Logger class.
- `PR #58 <https://github.com/martimunicoy/peleffy/pull/58>`_: Adds tests to validate consistency between PDB residue name and molecule tag.
- `PR #59 <https://github.com/martimunicoy/peleffy/pull/59>`_: Adds tests for the new conformer setter.


0.3.0 - Rotamers, OPLS2005, SMILES and stability improvements
-------------------------------------------------------------

This is a minor release that includes a refactoring of the classes and methods that involve the rotamer library builder. Besides, now it is possible to combine parameters from OPLS2005 and OFF. This release also contains a new method to define a molecule through a SMILES tag. It is still a preliminary version of the Open Force Field to PELE package which is under development.

New features
""""""""""""
- `PR #28 <https://github.com/martimunicoy/peleffy/pull/28>`_: Adds a new method to define a `Molecule` object through a SMILES tag. This molecule can be written as a PDB file later for PELE.
- `PR #31 <https://github.com/martimunicoy/peleffy/pull/31>`_: Adds the possibility to combine nonbonding and solvent parameters from OPLS2005 with bonding parameters from OFF.
- `PR #36 <https://github.com/martimunicoy/peleffy/pull/36>`_: Minor changes to improve the quality of the code.
- `PR #38 <https://github.com/martimunicoy/peleffy/pull/38>`_: Adds a new partial charge calculator that uses OPLS2005 to assign partial charges. Includes new flags in the CLI from main.py to combine bonding and nonbonding parameters and partial charges from OPLS2005.
- `PR #42 <https://github.com/martimunicoy/peleffy/pull/42>`_: Improves the documentation, adding a section specific for CLI-usage and API examples.
- `PR #46 <https://github.com/martimunicoy/peleffy/pull/46>`_: Adds a tag to Molecule class. Besides, the handling of Molecule names is improved. Both attributes can be set when initiating the molecule.

Bugfixes
""""""""
- `PR #22 <https://github.com/martimunicoy/peleffy/pull/22>`_: Fixes many bugs. For example, the default output name of the solvent parameters template is changed to `ligandParams.txt`, which is the name that PELE expects.
- `PR #32 <https://github.com/martimunicoy/peleffy/pull/32>`_: Minor fixes in ToolkitWrapper classes.
- `PR #34 <https://github.com/martimunicoy/peleffy/pull/34>`_: Improves the translation of dihedrals coming from the Open Force Fielf Toolkit and corrects the lack of exclusions in PELE 1-4 list that result from Impact's dihedral definitions.
- `PR #46 <https://github.com/martimunicoy/peleffy/pull/46>`_: Prevents molecule to be untagged when loading it from a SMILES tag.

Tests added
"""""""""""
- `PR #31 <https://github.com/martimunicoy/peleffy/pull/31>`_: Adds tests to validate some functions of the new SchrodingerToolkitWrapper.
- `PR #34 <https://github.com/martimunicoy/peleffy/pull/34>`_: Adds tests to further validate the assignment of parameters from the Open Force Field Toolkit.
- `PR #38 <https://github.com/martimunicoy/peleffy/pull/38>`_: Adds tests to validate the new OPLS charge calculator.
- `PR #46 <https://github.com/martimunicoy/peleffy/pull/46>`_: Adds tests to validate the name and tag assignment to Molecule class.


0.2.1
-----

This is a micro release that includes new features and parameters to configurate the behaviour of the program.
It is designed to be employed to run the first benchmarks of the implementation in PELE.
It also includes many stability improvements and an extended test coverage.

New features
""""""""""""
- `PR #15 <https://github.com/martimunicoy/peleffy/pull/15>`_: Adds a new method (Antechamber's gasteiger) to calculate partial charges.
- `PR #19 <https://github.com/martimunicoy/peleffy/pull/19>`_: Adds a new option to ignore terminal rotatable bonds of each rotamer's branch.
- `PR #17 <https://github.com/martimunicoy/peleffy/pull/17>`_: Adds and updates the documentation. However, it is still not completed.

Bugfixes
""""""""
- `PR #18 <https://github.com/martimunicoy/peleffy/pull/18>`_: Fixes some problems with proper and improper constructors.

Tests added
"""""""""""
- `PR #15 <https://github.com/martimunicoy/peleffy/pull/15>`_: Adds tests ensuring that the run_peleffy call from main and the partial charge calculators work as expected.
- `PR #19 <https://github.com/martimunicoy/peleffy/pull/19>`_: Adds tests to validate the construction of the `RotamerLibrary` class and the filtering of terminal rotatable bonds.


0.2.0
-----

This is a preliminary version of the Open Force Field to PELE package.

New features
""""""""""""

A first implementation of the package that allows to:

- Build a rotamer library for a small molecule using RDKit's API
- Build a template with the Molecular Mechanics' parameters for a small molecule using the Open Force Field Toolkit
- Assign the OBC implicit solvent parameters to a small molecule using the Open Force Field Toolkit
