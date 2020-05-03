.. _installationguide:

Installation guide
==================

C++
---

dCGP is a header-only library which has the following third party dependencies:

* `Boost <http://www.boost.org/>`_, various C++ utilities.
* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_, linear algebra library.
* `Pagmo <https://github.com/esa/pagmo2>`_, parallel optimization library.
* `tbb <https://github.com/intel/tbb>`_, lets you easily write parallel C++ programs that take full advantage of multicore performance.
* `AuDi <http://darioizzo.github.io/audi/>`_, high order automated differentiation library.
* `Symengine <https://github.com/symengine/symengine>`_, symbolic manipulation of math expressions.

After making sure the dependencies above are installed and found in your system, you may download
the latest dCGP code via git:

.. code-block:: console

   $ git clone https://github.com/darioizzo/dcgp.git

and configure your build using CMake. When done, type (in your build directory):

.. code-block:: console

   $ make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. 
To check that all went well compile the :ref:`quick-start example <getting_started_c++>`.

-----------------------------------------------------------------------

Python
------
The main functionalities of dCGP are exposed into a Python module called ``dcgpy`` which
can be installed from conda (OSx, linux and Win), pip (only linux) or by building the module.

Installing with conda
^^^^^^^^^^^^^^^^^^^^^
``dcgpy`` is available in the `conda <https://conda.io/en/latest/>`__ package manager
from the `conda-forge <https://conda-forge.org/>`__ channel. A single package is available:

* `dcgp-python <https://anaconda.org/conda-forge/dcccgp-python>`__, which contains the ``dcgpy`` python module.

In order to install ``dcgpy`` via conda, you just need
to add ``conda-forge`` to the channels:

.. code-block:: console

   $ conda config --add channels conda-forge
   $ conda install dcgp-python

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions
on how to setup and manage your conda installation.

Installing with pip (deprecated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We also provide the pip packages (mainly for linux 64 bit architectures and older versions).
Check on the `PyPi dcgpy page <https://pypi.org/project/dcgpy/>`_ if the needed package is provided.

.. code-block:: console

   $ pip install dcgpy

Building
^^^^^^^^^^^^^^^^^^^^^^^^^^

To build the module you need to have the Boost Python libraries installed and to activate the BUILD_DCGPY option from within CMake (and deselect BUILD_DCGP)

Check carefully what Python version is detected and what libraries are linked to. In particular select the correct boost_python
according to the Python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and Python module after install.

When done, type (in your build directory):

.. code-block:: console

   $ make install

To check that all went well fire-up your Python console and try the example in :ref:`quick-start example <getting_started_py>`.
