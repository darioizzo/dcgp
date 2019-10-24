Find an exact model for the Koza quintic problem
==================================================

In this first tutorial we show how to find an exact formula for some input data that do not require
any real valued constant. This is the easiest case for a symbolic regression task and thus makes it for a perfect entry tutorial.

We use the classic problem Koza quintic polynomial, that is x - 2x^3 + x^5.

Code:
^^^^^^^^
.. literalinclude:: ../../../examples/symbolic_regression_1.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output will be different on your computers as its non deterministic.

.. code-block:: python

    Gen:        Fevals:          Best:  Constants:    Formula:
       0              0        3898.35           []    [2*x0**3] ...
     500           2000        638.426           []    [x0**5] ...
    1000           4000        138.482           []    [(-x0**2 + x0**4)*x0] ...
    1500           6000        101.734           []    [-x0 + (-x0**2 + x0**4)*x0] ...
    1698           6792     5.2071e-30           []    [x0*(1 - x0**2) - x0**3*(1 - x0**2)] ...
    Exit condition -- ftol < 1e-08
    
    Best fitness: [5.2071e-30]
    Chromosome: [2, 0, 0, 3, 1, ... ]
    Pretty Formula: [(((((x0*x0)/(x0*x0))-(x0*x0))*x0)-(((((x0*x0)/(x0*x0))-(x0*x0))*x0)*(x0*x0)))]
    Prettier Formula: x0*(1 - x0**2) - x0**3*(1 - x0**2)
    Expanded Formula: x0 - 2*x0**3 + x0**5