dcgp::expression_weighted, A weighted Cartesian Genetic Program (dCGP-W)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class represents a **Weighted Cartesian Genetic Program**. Each node connection is associated to a weight so that more generic mathematical expressions
can be represented. When instantiated with the type *gdual<T>*, also the weights are defined as gduals, hence the program output can be expanded also with respect to the weights
thus allowing to train the weights using algorithms such as stochastic gradient descent, while the rest of the expression remains fixed. 


The class template can be instantiated using the types *double* or *gdual<T>*. 

.. figure:: ../_static/expression_weighted.png
   :alt: weighted dCGP expression
   :align: center

   A weighted dCGP expression

.. doxygenclass:: dcgp::expression_weighted
   :project: dCGP
   :members: