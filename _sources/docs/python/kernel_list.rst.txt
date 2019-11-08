Available kernels
----------------------------------

When constructing a :class:``dcgpy.kernel_set_double``, :class:``dcgpy.kernel_set_gdual_double `` 
or :class:``dcgpy.kernel_set_gdual_vdouble`` we can use the following names to add the corresponding
kernels to the set.

---------------------------------------------------------------------------

.. cssclass:: table-bordered table-striped

   +----------------+-----------------------+
   |Kernel name     |  Function             |
   +================+=======================+
   |"sum"           |addition               |
   +----------------+-----------------------+
   |"diff"          |subtraction            |
   +----------------+-----------------------+
   |"mul"           |multiplication         |
   +----------------+-----------------------+
   |"div"           |division               |
   +----------------+-----------------------+
   |"pdiv"          |protected division     |
   +----------------+-----------------------+
   |"sin"           |sine                   |
   +----------------+-----------------------+
   |"cos"           |cosine                 |
   +----------------+-----------------------+
   |"log"           | natural logarithm     |
   +----------------+-----------------------+
   |"exp"           |exponential            |
   +----------------+-----------------------+
   |"gaussian"      |gaussian               |
   +----------------+-----------------------+
   |"sqrt"          |square root            |
   +----------------+-----------------------+
   | **Suitable for dCGPANN**               |
   +----------------+-----------------------+
   |"sig"           |sigmoid                |
   +----------------+-----------------------+
   |"tanh"          |hyperbolic tangent     |
   +----------------+-----------------------+
   |"ReLu"          |rectified linear unit  |
   +----------------+-----------------------+
   |"ELU"           |exp linear unit        |
   +----------------+-----------------------+
   |"ISRU"          |sigmoid                |
   +----------------+-----------------------+
   |"sum"           |addition               |
   +----------------+-----------------------+
