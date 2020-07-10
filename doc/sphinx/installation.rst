.. _installationguide:

Installation guide
==================

C++
---

dCGP is a header-only library which has the following third party dependencies:

* `Boost <http://www.boost.org/>`_, various C++ utilities. (>=1.72).
* `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_, linear algebra library. (>=3.3.0)
* `Pagmo <https://github.com/esa/pagmo2>`_, parallel optimization library. (>=2.15).
* `tbb <https://github.com/intel/tbb>`_, lets you easily write parallel C++ programs that take full advantage of multicore performance. (>=2020.1).
* `AuDi <http://darioizzo.github.io/audi/>`_, high order automated differentiation library. (>=1.8).
* `obake <https://github.com/bluescarni/obake>`,  symbolic manipulation of sparse polynomials. (>=0.6.0). 
* `Symengine <https://github.com/symengine/symengine>`_, symbolic manipulation of math expressions. (>=0.6).

In case you are familiar with the conda package manager an environent ready for the dcgp installation can be created via the single command

.. code-block:: console

   $ conda create -n build_dcgp cmake cxx-compiler boost-cpp pagmo-devel tbb-devel audi symengine obake-devel

After making sure the dependencies above are installed and found in your system, you may download
the latest dCGP code via git:

.. code-block:: console

   $ git clone https://github.com/darioizzo/dcgp.git

and configure your build using CMake. For example, if you are using conda, make a build directory and type:

.. code-block:: console

   $ cmake -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_INSTALL_PREFIX=~/miniconda3/envs/build_dcgp/ -DCMAKE_PREFIX_PATH=~/miniconda3/envs/build_dcgp/ ../

where the conda environment has been assumed in ```~/miniconda3/envs/build_dcgp/```.

When done, type (in your build directory):

.. code-block:: console

   $ make 

When finished, to run the tests type:

.. code-block:: console

   $ make test

If succesfull, you may now install cgp:

.. code-block:: console

   $ make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. 
To check that all went well compile the :ref:`quick-start example <getting_started_c++>`.



-----------------------------------------------------------------------

Python
------
The main functionalities of dCGP are exposed into a Python module called ``dcgpy`` which
can be installed from conda (OSx, linux and Win), pip (only linux) or by building the module.

The following third party dependencies are required to have full access to the ``dcgpy`` API:

* `numpy <https://numpy.org/>`_, The fundamental package for scientific computing with Python. (>=1.18)
* `matplotlib <https://matplotlib.org/>`_,  A comprehensive library for creating static, animated, and interactive visualizations in Python. (>=3.2)
* `pyaudi <http://darioizzo.github.io/audi/>`_, A library that implements the differential algebra of Taylor truncated polynomials. (>=1.8)
* `sympy <https://www.sympy.org/en/index.html>`_, A Python library for symbolic mathematics. (>=1.6)
* `graphviz <https://graphviz.readthedocs.io/en/stable/>`_, A simple pure-Python interface for the Graphviz graph-drawing software. (>=2.42)


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

note that all the required dependencies will be installed automatically, as well as the C++ ```dcgp``` headers.

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions
on how to setup and manage your conda installation.

You may test the successfull installation by running the python tests typing:

.. code-block:: console

   $ python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"


Building
^^^^^^^^^^^^^^^^^^^^^^^^^^

To build the python module you need to first install the dcgp C++ header library and its dependencies (see above) as well as the additional dependency:

* `pybind11 <https://github.com/pybind/pybind11>`, Seamless operability between C++11 and Python. (>=2.5.0). 

In case you are familiar with the conda package manager an environent ready for the python module installation can be created via the single command

.. code-block:: console

   $ conda create -n build_dcgp cmake cxx-compiler boost-cpp pagmo-devel tbb-devel audi symengine obake-devel pybond11

Install the latest dCGP code via git:

.. code-block:: console

   $ git clone https://github.com/darioizzo/dcgp.git

After installing the C++ dcgp library (see above) and making sure your environment is correctly set up to find all dependencies, you can 
configure your python module build using CMake. For example, if you are using conda, in the build directory type:

.. code-block:: console

   $ cmake -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_INSTALL_PREFIX=~/miniconda3/envs/build_dcgp/ -DCMAKE_PREFIX_PATH=~/miniconda3/envs/build_dcgp/ -DDCGP_BUILD_DCGP=OFF -DDCGP_BUILD_DCGPY=ON ../

where the conda environment has been assumed in ```~/miniconda3/envs/build_dcgp/```.

When done, type (in your build directory):

.. code-block:: console

   $ make install

To check that all went well fire-up your Python console and try the example in :ref:`quick-start example <getting_started_py>`.
