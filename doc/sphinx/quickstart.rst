.. quickstart examples


Quick start examples
====================

.. contents::


C++
---

After following the :ref:`installationguide` you will be able to compile and run your first C++ dCGP program:

.. _getting_started_c++:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

Place it into a getting_started.cpp text file and compile it with:

.. code-block:: bash

   g++ -std=c++11 getting_started.cpp -lmpfr -lgmp -pthread

-----------------------------------------------------------------------

Python
------

If you have successfully compiled and installed dcgpy following the :ref:`installationguide` you will be able to test its use by running the following script:

.. _getting_started_py:

.. literalinclude:: ../../doc/examples/getting_started.py
   :language: python
   :linenos:

Place it into a getting_started.py text file and run it with:

.. code-block:: bash

   python getting_started.py

We recommend the use of Jupyter or Ipython do enjoy dcgpy the most.

