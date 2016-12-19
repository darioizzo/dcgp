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

.. _notebooks:

Notebooks
^^^^^^^^^

Follow the links below to visualize Jupyter notebooks on the use of dCGP:

- `Learning constants in a symbolic regression task I <https://github.com/darioizzo/d-CGP/blob/master/examples/learning_constants.ipynb>`_ (by Dario Izzo)

- `Learning constants in a symbolic regression task II <https://github.com/darioizzo/d-CGP/blob/master/examples/learning_constants2.ipynb>`_ (by Dario Izzo)

- `Weighted dCGP for a symbolic regression task <https://github.com/darioizzo/d-CGP/blob/master/examples/weighted_symbolic_regression.ipynb>`_ (by Alessio Mereta)

- `Solving differential equations with dCGP <https://github.com/darioizzo/d-CGP/blob/master/examples/solving_odes.ipynb>`_ (by Dario Izzo)

- `Discovery of prime integrals <https://github.com/darioizzo/d-CGP/blob/master/examples/finding_prime_integrals.ipynb>`_ (by Dario Izzo)
