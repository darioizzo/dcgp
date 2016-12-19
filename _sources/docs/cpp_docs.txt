.. cpp docs

C++ Documentation
=================

.. contents::

Classes
-------

expression: a Cartesian Genetic Program
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/expression.png
   :alt: dCGP expression
   :align: center

   A dCGP expression

.. doxygenclass:: dcgp::expression
   :project: dCGP
   :members:

----------------------------------------------------------

expression_weighted: a weighted Cartesian Genetic Program
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: ../_static/expression_weighted.png
   :alt: weighted dCGP expression
   :align: center

   A weighted dCGP expression

.. doxygenclass:: dcgp::expression_weighted
   :project: dCGP
   :members:

----------------------------------------------------------

kernel: a function defining the generic CGP node
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: dcgp::kernel
   :project: dCGP
   :members:

----------------------------------------------------------

kernel_set: a set of kernels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _kernels:

.. cssclass:: table-bordered table-striped

   +----------------+-----------------------+
   |Kernel name     |  Function             |
   +================+=======================+
   |"sum"           |addition               |
   +----------------+-----------------------+
   |"diff"          |subtraction            |
   +----------------+-----------------------+
   |"mul"           |multiplication         |
   +----------------+-----------------------+
   |"div"           |division               |
   +----------------+-----------------------+
   |"pdiv"          |protected division     |
   +----------------+-----------------------+
   |"sig"           |sigmoid                |
   +----------------+-----------------------+
   |"sin"           |sine                   |
   +----------------+-----------------------+
   |"cos"           |cosine                 |
   +----------------+-----------------------+
   |"log"           |logarithm              |
   +----------------+-----------------------+
   |"exp"           |exponential            |
   +----------------+-----------------------+

.. doxygenclass:: dcgp::kernel_set
   :project: dCGP
   :members:
