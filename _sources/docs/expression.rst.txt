dcgp::expression, The standard Cartesian Genetic Program (dCGP)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This class represents a **Cartesian Genetic Program**. Since that is, essentially, an artificial genetic encoding for a mathematical expression, we named the templated class *expression*.
The class template can be instantiated using the types *double* or *gdual<T>*. In the case of *double*, the class would basically reproduce a canonical CGP expression. In the case of *gdual<T>*
the class would operate in the differential algebra of truncated Taylor polynomials with coefficients in *T*, and thus provide also any order derivative information on the program 
(i.e. the Taylor expansion of the program output with respect to its inputs).

.. figure:: ../_static/expression.png
   :alt: dCGP expression
   :align: center

   A dCGP expression

.. doxygenclass:: dcgp::expression
   :project: dCGP
   :members: