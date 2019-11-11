Find an exact model inclduing one parameter using a memetic approach
====================================================================

In this second tutorial we show how to find a model for our input data when we also
want to learn some constants. 

Constants can, in general, be learned via two main techniques:
 1. evolutionary (common and standard practice in GP)
 2. memetic (original with dCGP)

The difference is that the evolutionary approach cannot find the correct and exact values for constants, only approximations. 
In this tutorial we follow the evolutionary approach 2. In the next tutorial we will follow a memetic approach 2.
We use the problem P1 from the dcgp::gym, that is x^5 - pi*x^3 + x

Code:
^^^^^^^^
.. literalinclude:: ../../../../examples/symbolic_regression_3.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output is non deterministic. Sometimes, with a bit of luck :)
the problem is solved exaclty (loss goes to zero). The following output reports one of these occurences.
Note that in the traditional evolutionary approach this result is incredibly hard to obtain.

.. code-block:: python

   Gen:        Fevals:          Best:	Constants:	Formula:
      0              0        4009.59	[-2.41031]	[(c1 + x0)*x0] ...
     50            200       0.978909	[0.238557]	[x0**2*c1*(-x0 + x0**4) - (c1 + x0 + x0* ...
    100            400        0.84565	[0.240548]	[x0**2*c1*(-x0 + x0**4) - (c1 + 2*x0)] ...
    150            600       0.761757	[0.240032]	[-2*x0 + x0**2*c1*(-x0 + x0**4)] ...
    200            800      0.0170582	[-1.16484]	[(-x0 + x0**3)*(c1 + x0**2) - x0**3] ...
    250           1000      0.0170582	[-0.164837]	[(-x0 + x0**3)*(c1 + x0**2) - (-x0 + 2*x ...
    300           1200      0.0170582	[-0.164837]	[(-x0 + x0**3)*(c1 + x0**2) - (-x0 + 2*x ...
    350           1400      0.0170582	[-0.164837]	[(-x0 + x0**3)*(c1 + x0**2) - (-x0 + 2*x ...
    357           1428    2.17578e-29	[-1.14159]	[x0**3*(c1 + x0**2) - (-x0 + 2*x0**3)] ...
   Exit condition -- ftol < 1e-08
   
   Best fitness: [2.17578e-29]
   Chromosome: [-1.14159, 2, 0, 0, 2, ... ]
   Pretty Formula: [(((c1+(x0*x0))*((x0*x0)*x0))-(((x0*x0)*x0)+(((x0*x0)*x0)-x0)))]
   Prettier Formula: x0**3*(c1 + x0**2) - (-x0 + 2*x0**3)
   Expanded Formula: x0 + x0**3*c1 - 2*x0**3 + x0**5
