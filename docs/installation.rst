.. _installation ::

Installation
************

Installing via `conda`
======================
The more straightforward way to install `peleffy` along with the required
dependencies is through the `conda <http://www.continuum.io/blog/conda>`_
package manager.

First, you need to install the
`miniconda <http://conda.pydata.org/miniconda.html>`_ distribution, which is
the minimal installation of the Anaconda Python package.

To install the Python 3 version on ``linux`` (on ``bash`` systems):

.. code-block:: bash

   $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   $ bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

It can also be installed on ``osx`` with:

.. code-block:: bash

   $ curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O
   $ bash ./Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda3

Once installed, the bases environment of `conda` can be loaded with
the following commands:

.. code-block:: bash

   $ source ~/miniconda3/etc/profile.d/conda.sh
   $ conda activate base

You can also create a custom conda environment to handle peleffy:

.. code-block:: bash

   $ conda create --name peleffy_env

which can be activated and deactivated with the two commands from below:

.. code-block:: bash

   $ conda activate peleffy_env
   $ conda deactivate

To install de dependencies of peleffy, the following `conda` channels need
to be added and updated:

.. code-block:: bash

   $ conda config --add channels omnia --add channels conda-forge --add channels martimunicoy
   $ conda update --all

Finally, you can install the latest stable build of `peleffy`

.. code-block:: bash

   $ conda install peleffy


Installing via `PyPI`
=====================

`peleffy` can also be installed through `PyPI <https://pypi.org>`_
with the following command:

.. code-block:: bash

   $ pip install peleffy

However, with `PyPI` some of the required dependencies of `peleffy` are not
installed and have to be installed manually such as:

- Open Force Field Toolkit
- RDKit
- AmberTools

For this reason, the installation through `conda` is recommended.


External dependencies
=====================

Some of the functionalities of `peleffy` require external dependencies.
They are normally included with the standard `conda` installation, as
explained above. However, the Schrodinger toolkit must be installed
manually. It is only required when combining `Open Force Field` parameters
with `OPLS2005` (as it uses the Schrodinger's `ffld_server`). Nevertheless,
in case that Schrodinger dependencies are missing, `peleffy` can still be
employed to generate pure `Open Force Field` parameters.

The easiest way to get a valid Schrodinger installation is downloading
`Free Maestro <https://www.schrodinger.com/freemaestro>`_. It can be
installed in both platforms that are supported by `peleffy`: Linux and
MacOS. Once installed, `peleffy` will need an environment variable to be
set in order to known the Schrodinger's installation path. So, please,
check that the following environment variable is set before running
`peleffy` if you plant to work with `OPLS2005` parameters:

.. code-block:: bash

   $ export SCHRODINGER=/path/to/Schrodinger/installation/

For example, in MacOS, a typical installation path is
`/opt/schrodinger/suites2020-2/`. Therefore:

.. code-block:: bash

   $ export SCHRODINGER=/opt/schrodinger/suites2020-2/

This variable must be set every time `peleffy` is employed to work with
`OPLS2005` parameters in a new console session. To avoid future concerns about
this issue, you can set the environment variable automatically every time you
initiate a bash session in your console. You can do so by modifying your
`.bashrc`, `.bash_profile` or `.zshrc` (in case of a `zsh` shell) by adding the
line above.
