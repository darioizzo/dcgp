.. python docs

Python Documentation
====================

.. contents::


Types of dCGPs
--------------

expression_double
^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_double
    :members:

expression_gdual_double
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_gdual_double
    :members:

expression_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_gdual_vdouble
    :members:

expression_weighted_double
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_weighted_double
    :members:

expression_weighted_gdual_double
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_weighted_gdual_double
    :members:

expression_weighted_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.expression_weighted_gdual_vdouble
    :members:

expression_ann_double
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. important::
   This Cartesian Genetic Program is able to encode an Artificial Neural Network. Weights and biases are added to the acyclic graph
   as well as extra methods to allow to perform backpropagation (in parallel), to visualize the network and more ...

.. autoclass:: dcgpy.expression_ann_double
    :members:

Non linearities
--------------------

kernel_double
^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_double
    :members:

kernel_gdual_double
^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_gdual_double
    :members:

kernel_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_gdual_vdouble
    :members:

kernel_set
^^^^^^^^^^

For a list of the available kernels see :ref:`kernels <kernels>`.

kernel_set_double
^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_set_double

    .. automethod:: dcgpy.kernel_set_double.push_back()

kernel_set_gdual_double
^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_set_gdual_double

    .. automethod:: dcgpy.kernel_set_gdual_double.push_back()

kernel_set_gdual_vdouble
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: dcgpy.kernel_set_gdual_vdouble

    .. automethod:: dcgpy.kernel_set_gdual_vdouble.push_back()
