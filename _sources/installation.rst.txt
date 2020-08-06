.. _installation ::

Installation
************

Installing via `conda`
======================
The more straightforward way to install `offpele` along with the required
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

You can also create a custom conda environment to handle offpele:

.. code-block:: bash

   $ conda create --name offpele_env

which can be activated and deactivated with the two commands from below:

.. code-block:: bash

   $ conda activate offpele_env
   $ conda deactivate

To install de dependencies of offpele, the following `conda` channels need
to be added and updated:

.. code-block:: bash

   $ conda config --add channels omnia --add channels conda-forge --add channels martimunicoy
   $ conda update --all

Finally, you can install the latest stable build of `offpele`

.. code-block:: bash

   $ conda install offpele


Installing via `PyPI`
=====================

`offpele` can also be installed through `PyPI <https://pypi.org>`_
with the following command:

.. code-block:: bash

   $ pip install offpele

However, with `PyPI` some of the required dependencies of `offpele` are not
installed and have to be installed manually such as:

- Open Force Field Toolkit
- RDKit
- AmberTools

For this reason, the installation through `conda` is recommended.
