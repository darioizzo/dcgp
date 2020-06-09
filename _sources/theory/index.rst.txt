Background
====================

Cartesian Genetic Programming (CGP) is a tree based representation of a computer program introduced by 
Julian F. Miller and Peter Thomson in 1997. Several open source codes exist that implement some
form of CGP, Miller himself maintain a resource page `here <https://www.cartesiangp.com/resources>`_.

The dcgpy library here introduced, differs in several aspects from all existing implementations:
it allows to perform any-order derivatives on the represented programs (hence implementing a
differentiable form of genetic programming), it allows to use Python to use runtime scripting, 
it is thread-safe and it allows to define kernels as generic functors.

.. toctree::
  :maxdepth: 1

  cgp
  dcgp
