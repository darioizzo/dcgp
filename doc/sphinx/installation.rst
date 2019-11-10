.. _installationguide:

Installation guide
==================

C++
---

dCGP is a header-only library which has the following third party dependencies:

* `AuDi <http://darioizzo.github.io/audi/>`_, header-only library (git clone https://github.com/darioizzo/audi.git)
* `Boost <http://www.boost.org/>`_, headers only (needs the libraries if you intend to use the python bindings)

After making sure the dependencies above are installed in your system, you may download the latest dCGP version via git:

.. code-block:: bash

   git clone https://github.com/darioizzo/d-CGP.git

and configure your build using CMake. When done, type (in your build directory):

.. code-block:: bash

   make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. To check that all went well compile the :ref:`quick-start example <getting_started_c++>`.

-----------------------------------------------------------------------

Python
------
The main functionalities of dCGP are exposed into a Python module called dcgpy.
It can be installed either directly from pip or by building the module.

Installing with pip
^^^^^^^^^^^^^^^^^^^
On a Win 64bit system or a Linux based system (32 or 64 bits), the Python package dcgpy (Python binding of the C++ code) can be installed using ``pip``:

.. code-block:: bash

   pip install dcgpy

Building the python module
^^^^^^^^^^^^^^^^^^^^^^^^^^

To build the module you need to have the Boost Python libraries installed and to activate the BUILD_DCGPY option from within CMake (and deselect BUILD_DCGP)

Check carefully what Python version is detected and what libraries are linked to. In particular select the correct boost_python
according to the Python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and Python module after install.

When done, type (in your build directory):

.. code-block:: bash

   make install

To check that all went well fire-up your Python console and try the example in :ref:`quick-start example <getting_started_py>`.
