.. _installationguide:


Installation guide
==================

.. contents::


C++
---

AuDi is a header only library which has the following third party dependencies

* The `boost <http://www.boost.org/>`_ C++ libraries: the following boost libraries are necessary: boost_system, boost_unit_test_framework, boost_timer, boost_chrono. boost headers must be found in the system
* `piranha <http://bluescarni.github.io/piranha/index.html>`_: piranha headers must be found in the system
* `GNU MPFR library <http://www.mpfr.org/>`_: needed by piranha
* `The GMP multiprecision library <https://gmplib.org/>`_: needed by piranha

After making sure the dependencies above are installed in your system (most linux / osx package managers include them), you may download the latest AuDi version via git:

.. code-block:: bash

   git clone https://github.com/darioizzo/audi.git

and onfigure your build using CMake. When done, type (in your build directory):

.. code-block:: bash

   make install

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. To check that all went well compile the :ref:`quick-start example <getting_started>`.

-----------------------------------------------------------------------

Python
------

The main functionalities of AuDi are exposed into a python module called pyaudi. To create the module you need to have
the boost python libraries installed and activate the BUILD_PYAUDI option from within cmake.

Check carefully what python version is detected and what libraries are linked to. In particular select the correct boost_python
according to the python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and python module after install.

When done, type (in your build directory):

.. code-block:: bash

   make install

To check that all went well fire-up your python console and try the example in :ref:`quick-start example <getting_started>`.
