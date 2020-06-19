.. quickstart examples


Quick start examples
====================

C++
---

After following the :ref:`installationguide` you will be able to compile and run your first C++ dCGP program, 
put the following text into a ``getting_started.cpp`` file:

.. _getting_started_c++:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

To compile it, create also, in the same folder, a ``CmakeLists.txt`` file with the following content:

.. code-block:: cmake

    project(getting_started)

    cmake_minimum_required(VERSION 3.2)
    cmake_policy(SET CMP0057 NEW)

    find_package(dcgp REQUIRED)

    add_executable(getting_started getting_started.cpp)
    target_link_libraries(getting_started Dcgp::dcgp)

    set_property(TARGET getting_started PROPERTY CXX_STANDARD 17)
    set_property(TARGET getting_started PROPERTY CXX_STANDARD_REQUIRED YES)
    set_property(TARGET getting_started PROPERTY CXX_EXTENSIONS NO)

then:

.. code-block:: console

    $ mkdir build
    $ cd build
    $ cmake ../
    $ make
    $ ./getting_started

-----------------------------------------------------------------------

Python
------

If you have successfully compiled and installed dcgpy following the :ref:`installationguide` you will be able to test its
use by running the following script:

.. _getting_started_py:

.. literalinclude:: ../../doc/examples/getting_started.py
   :language: python
   :linenos:

Place it into a ``getting_started.py`` text file and run it with:

.. code-block:: console

   $ python getting_started.py

We recommend the use of `Jupyter <https://jupyter.org/>`_ or `Ipython <https://github.com/ipython/ipython>`_ to enjoy ``dcgpy`` the most.

