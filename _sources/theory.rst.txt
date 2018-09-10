.. contents::

Cartesian Genetic Programming
=============================
Cartesian Genetic Programming (CGP) is a form of Genetic Programming where the program
representation consists of a graph. For more information visit
`Cartesian Genetic Programming <http://www.cartesiangp.co.uk/>`_.

dCGP
====
Differentiable Cartesian Genetic Programming (dCGP) is a recent development of CGP
that adds the information about the derivatives of the output nodes (the programs,
or expressions encoded) with respect to the input nodes (the input values) and/or the
weights (that can be associated with every connection in the graph, similarly to Neural Networks).

.. figure:: _static/expression_theory.png
   :alt: weighted dCGP expression
   :align: center
   :width: 800px

   A weighted dCGP expression

The derivatives are obtained by means of Automatic Differentiation through the
`AuDi <http://darioizzo.github.io/audi/>`_ library, that implements the Taylor truncated
polynomial algebra. See :ref:`reference` for more details.

The use of the derivatives of the outputs (and hence of any fitness function that is a
combination of these) with respect to inputs/weights enables a number of new applications
currently the subject of active research.

The evolution of the genetic program can now be supported by using the information
on the derivatives, hence enabling for the equivalent of back-propagation in Neural Networks.

Furthermore, the fitness function itself can be defined in terms of the derivatives,
allowing for additional tasks beyond simple regression, *e.g.*:

* solving differential equations,
* learning differential models,
* capturing conserved quantities in dynamical systems.

.. _reference:

Reference
^^^^^^^^^^

Dario Izzo, Francesco Biscani, and Alessio Mereta. `Differentiable Genetic Programming. <https://arxiv.org/pdf/1611.04766v1.pdf>`_ arXiv preprint arXiv:1611.04766 (2016).
