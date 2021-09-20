PELE Force Field Yielder
========================

The `peleffy` (PELE Force Field Yielder) is a Python package that builds
PELE-compatible force field templates. The current supported force fields
are:

* The following force field from the `Open Force Field toolkit <https://github.com/openforcefield/openforcefield>`_:

  * openff_unconstrained-2.0.0.offxml
  
  * openff_unconstrained-1.3.0.offxml

  * openff_unconstrained-1.2.1.offxml

  * openff_unconstrained-1.2.0.offxml

  * openff_unconstrained-1.1.1.offxml

  * openff_unconstrained-1.1.0.offxml

  * openff_unconstrained-1.0.1.offxml

  * openff_unconstrained-1.0.0.offxml

* OPLS2005

It also generates other files required by PELE such as rotamer libraries or
solvent templates.


User guide
----------

.. toctree::
  :maxdepth: 1

  installation
  releasehistory
  usage
  examples


API documentation
-----------------

.. toctree::
  :maxdepth: 1

  topology
  forcefield
  template
  solvent
