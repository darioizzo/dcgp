.. cpp docs

C++ Documentation
=================

Types of dCGPs
--------------------

Several types of Cartesian Genetic Program are provided in dCGP. Since the inner computational graph present in each CGP defines,
some kind of mathematical expression, we use the term expression in all names for the different classes. 
We provide three types of CGPs:

* **expression**: this is the original CGP as introduced by Miller in 1999.
* **expression_weighted**: this adds to the original CGP formulation weights on each of the graph edges. (original with dCGP - 2016)
* **expression_ann**: this represents an Artificial Neural Network inclusive of biases and weights, via a CGP and allows to learn the model parameters using backproagation. (original with dCGP - 2018)

Both **expression** and **expression_weighted** can operate over different numerical types and thus are templated. 
For example, **expression** can operate over double, in which case the result of evaluating the inner computational graph will be a float,
but also over generalized dual numbers (gduals), in which case the result of evaluating the inner computational graph will also be a gdual
(hence it will contain all the program derivatives up to the chosen order.)

Another important type some CGP can operate upon in the vectorized gdual. This type is the same as the gdual type, but its vectorized,
allowing order of magnitude speed ups when a CGP needs to be evaluated over several points (such as in the case aof a loss evaluation)

The type can be selcted via the template parameter

.. toctree::
  :maxdepth: 1

  expression
  expression_weighted
  expression_ann


Kernels
--------------------

Kernels (also called non-linearities in the ANN literature) describe the fundamental computational units of a CGP. 
Things like addition, multiplication, trigonomoetric functions are all kernels. The following classes allow the user
the definition of their own kernels able to operate on the choosen type. Note that some commonly used kernels are
provided. They can be used when constructing a *dcgp::kernel_set*.

.. toctree::
  :maxdepth: 1

  kernel
  kernel_set
  