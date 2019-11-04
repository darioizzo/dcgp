.. python docs

Python Documentation
====================

.. contents::


Types of dCGPs
--------------
Several types of Cartesian Genetic Program are provided in dcgpy. Since the inner computational graph defines, in all cases,
some kind of mathematical expression, we use the term expression in all names for the different classes. 
We essentially provide three types of dCGP:

* **expression**: this is the original CGP as introduced by Miller in 1999.
* **expression_weighted**: this adds to the original CGP formulation weights on each of the graph edges. (original with dCGP - 2016)
* **expression_ann**: this represents an Artificial Neural Network inclusive of biases and weights, via a CGP and allows to learn the model parameters using backproagation. (original with dCGP - 2018)

Each of the above CGPs can operate over different numerical types. For example **expression** can operate over floats, in which case the 
result of evaluating the inner computational graph will be a float, but also on gduals, in which case, the result of evaluating the
inner computational graph will be a gdual (hence it will contain all the program derivatives up to the chosen order.)

Another important type some CGP can operate upon is the vectorized gdual. This type is the same as the gdual type, but its vectorized,
allowing order of magnitude speed ups when a CGP needs to be evaluated over several points (such as in the case aof a loss evaluation)

.. toctree::
  :maxdepth: 1

  expression
  expression_weighted
  expression_ann

----------------------------------------------------------------------------------

kernels
---------

Kernels, (also called non-linearities in the ANN literature) describe the fundamental computational units of a CGP. 
Things like addition, multiplication, trigonomoetric functions are all kernels. The class ``kernel`` allow the user
the definition of their own kernels able to operate on the choosen type. Some popular kernels are shipped with dcgpy
and can be easily instantiated via the class ``kernel_set`` 

.. toctree::
  :maxdepth: 1

  kernel
  kernel_set

----------------------------------------------------------------------------------

 
Symbolic Regression Gym
-----------------------------

For the specific application of symbolic regression, we make available a number of data sets coming from 
various sources. These are collected in what its called the ``dcgpy gym`` and have all the same interface.

.. toctree::
  :maxdepth: 1

  koza_quintic
  p_problems
  vladislavleva_problems
  nist_problems

