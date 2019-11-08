.. cpp docs

C++ Documentation
=================

.. contents::

Kernels
--------- 

.. figure:: ../../_static/kernels.png
   :alt: dCGP kernels
   :align: right
   :scale: 40 %

**Kernels**, (also called **non-linearities** in the ANN literature) describe the fundamental computational units of a CGP. 
Things like addition, multiplication, trigonomoetric functions are all kernels. The templated class :cpp:class:`dcgp::kernel` allow 
the user the definition of their own kernels able to operate on the choosen type. 

The most popular kernels are already coded and shipped with dcgpy. A python list containing multiple kernels can be easily 
instantiated via the :cpp:class:`dcgp::kernel_set`.

.. toctree::
  :maxdepth: 1

  kernel
  kernel_set
  kernel_list

----------------------------------------------------------------------------------

Types of dCGPs
--------------

Several types of **Cartesian Genetic Program** are provided in dcgpy. Since a dCGP is some kind of a mathematical expression,
we use the term expression to name the related, different, classes. 

We essentially provide three types of CGPs:

* **expression**: this is the original CGP as introduced by Miller in 1999.
* **expression_weighted**: this adds to the original CGP formulation weights on each of the graph edges. (original with dCGP - 2016)
* **expression_ann**: this represents an Artificial Neural Network inclusive of biases and weights, via a CGP and allows to learn the model parameters using backproagation. (original with dCGP - 2018)

Each of the above CGPs can operate over different numerical types, hence the corresponding classes are templated. 
For example a :cpp:class:`dcgp::expression` can operate over doubles (``T`` = ``double``), in which case the 
result of evaluating the inner computational graph is a ``double``, but also on a ``gdual`` (``T`` = :cpp:class:`audi::gdual <gdual>` ``<Cf>``)
with ``Cf`` = ``double``, in which case, the result of evaluating the inner computational graph will
be a :cpp:class:`audi::gdual <gdual>` ``<Cf>`` and hence it will contain all the program derivatives
up to the chosen order and is thus referred to as a dCGP. 

Another important type some CGPs can operate upon is the :cpp:class:`audi::gdual <gdual>` ``<Cf>`` with ``Cf`` = ``vectorized_double``.
This type offers order of magnitude speed ups whenever a CGP needs derivatives and to be evaluated over several points 
(such as in the case of a loss evaluation).

.. toctree::
  :maxdepth: 1

  expression
  expression_weighted
  expression_ann

----------------------------------------------------------------------------------

 
Symbolic Regression
---------------------

.. figure:: ../../_static/non_linear_regression.png
   :alt: non-linear regression
   :align: right
   :scale: 70 %

Mathematically, a symbolic regression problem is a global optimization problem. In order to facilitate its solution,
a number of classes have been designed to interface a :cpp:class:`dcgp::expression` to the optimisation suite pygmo. 
In particular we provide UDPs (in pagmo's jargon user defined problems) that can be used to build
:cpp:class:`pagmo::problem <pagmo::problem>` objects (representing symbolic regression problems) and UDAs 
(in pagmo's jargon user defined algorithms) that can be used to build :cpp:class:`pagmo::algorithm <pagmo::algorithm>` objects.

.. toctree::
  :maxdepth: 1

  symbolic_regression
  es4cgp
  gd4cgp
  mes4cgp


We also make available, as a gym to test the capabilities of various proposed methodologies, a number of data sets
coming from different scientific sources. These are collected in what we call the ``dcgpy gym``.

.. toctree::
  :maxdepth: 1

  koza_quintic
  p_problems
  vladislavleva_problems
  nist_problems



  