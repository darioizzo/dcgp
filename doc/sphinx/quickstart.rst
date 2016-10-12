.. quickstart examples


Quick start examples
====================

.. contents::


C++
---

After following the :ref:`installationguide` you will be able to compile and run your first C++ AuDi program:

.. _getting_started:

.. literalinclude:: ../../doc/examples/getting_started.cpp
   :language: c++
   :linenos:

Place it into a getting_started.cpp text file and compile it with:

.. code-block:: bash

   g++ -std=c++11 getting_started.cpp -lmpfr -lgmp -pthread

-----------------------------------------------------------------------

Python
------

If you have succesfully compiled and installed pyaudi following the :ref:`installationguide` you will be able to test its use typing the following script.

.. literalinclude:: ../../doc/examples/getting_started.py
   :language: python
   :linenos:

Place it into a getting_started.py text file and run it with 

.. code-block:: bash

   python getting_started.py

We reccomend the use of Jupyter or ipython do enjoy pyaudi the most. 

Notebooks
^^^^^^^^^

Follow the links below to visualize juypiter notebooks on the use of pyaudi.

- `The very basics <https://github.com/darioizzo/audi/blob/master/examples/example00.ipynb>`_: (by Francesco Biscani and Dario Izzo)

- `Understanding gduals and floats <https://github.com/darioizzo/audi/blob/master/examples/example01.ipynb>`_: (Dario Izzo)

- `Training an artificial neural network <https://github.com/darioizzo/audi/blob/master/examples/example11.ipynb>`_: (by Carlos Sanchez)

- `Differential Intelligence <https://github.com/darioizzo/audi/blob/master/examples/example10.ipynb>`_: (by Dario Izzo)

