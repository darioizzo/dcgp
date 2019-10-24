Finding an entire non dominated front of formulas.
==================================================

In this fourth tutorial on symbolic regression we solve the multiobjective symbolic regression problem.
The Mean Squared Error (i.e the loss) of our model is considered next to the model complexity to determine how good 
a certain model is. The result is thus a whole non-dominated front of models.

This case is arguably the most complete and useful among  symbolic regression tasks. We use here
the problem vladi6 from the dcgp::gym, that is: 6*cos(x*sin(y))

Code:
^^^^^^^^
.. literalinclude:: ../../../examples/symbolic_regression_4.cpp
   :language: c++
   :linenos:

Output:
^^^^^^^
Note: the actual output will be different on your computers as its non deterministic.

.. code-block:: python

  Gen:        Fevals:     Best loss: Ndf size:  Compl.:
      0              0        4.10239        11        62
     10           1000        1.42703         7        74
     20           2000       0.828554         8        45
     30           3000       0.374803        13        78
     40           4000       0.164032        16        66
     50           5000        0.03552        14        48
     60           6000      0.0200792        11        45
     70           7000              0        10        19
     80           8000              0        10        19
     90           9000              0        10        19
    100          10000              0        11        19
    110          11000              0        10        18
    120          12000              0        11        18
    130          13000              0        11        18
    140          14000              0        11        18
    150          15000              0        11        18
    160          16000              0        11        18
    170          17000              0        11        18
    180          18000              0        11        18
    190          19000              0        11        18
    200          20000              0        11        18
    210          21000              0        11        18
    220          22000              0        11        18
    230          23000              0        11        18
    240          24000              0        10        18
    250          25000              0        10        18
    Exit condition -- max generations = 250
    
    Non dominated Front at the end:
    1  - Loss: 0            Complexity:    18   Formula:  c1*cos(x0*sin(x1))
    2  - Loss: 0.844272     Complexity:    17   Formula:  3*c1 - 2*x0 - 2*x0*x1
    3  - Loss: 1.10197      Complexity:    15   Formula:  4*c1 - 3*x0 - x0*x1
    4  - Loss: 1.17331      Complexity:    14   Formula:  3*c1 - 3*x0 - 2*x1
    5  - Loss: 1.27379      Complexity:    10   Formula:  c1 - 3*x0 - x1
    6  - Loss: 1.92403      Complexity:    9    Formula:  2*c1 - 3*x0
    7  - Loss: 1.92403      Complexity:    7    Formula:  c1 - 3*x0
    8  - Loss: 3.22752      Complexity:    5    Formula:  c1 - x0
    9  - Loss: 4.74875      Complexity:    2    Formula:  c1
    10 - Loss: 4.8741       Complexity:    1    Formula:  4
