dcgp::expression_ann, An artificial neural network Cartesian Genetic Program (dCGP-ANN)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class represents a **Artificial Neural Network Cartesian Genetic Program**. Each node connection is associated to a weight and each node to a bias. Only a subset of the kernel functions
is allowed, including the most used nonlinearities in ANN research: *tanh*, *sig*, *ReLu*, *ELU* and *ISRU*. The resulting expression can represent any feed forward neural network but also other
less obvious architectures. Weights and biases of the expression can be trained using the efficient backpropagation algorithm (gduals are not allowed for this class, they correspond to forward mode
automated differentiation which is super inefficient for deep networks ML.)

.. figure:: ../_static/expression_ann.png
   :alt: weighted dCGP expression
   :align: center

   A, small, artificial neural network as using the dCPP-ANN approach.

.. doxygenclass:: dcgp::expression_ann
   :project: dCGP
   :members:
