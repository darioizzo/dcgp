kernel_set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Classes such as :cpp:class:`dcgp::expression`, :cpp:class:`dcgp::symbolic_regression` and others are constructed 
from an ``std::vector`` of :cpp:class:`dcgp::kernel`. Assembling such a ``std::vector`` is an operation that is greatly 
facilitated by the class :cpp:class:`dcgp::kernel_set`.

Intended use of the class is:

.. highlight:: c++

.. code-block:: c++

   // (optional) - Say we have a custom kernel named f 
   kernel<double> f(my_sum<double>, print_my_sum, "my_sum");
   // We can construct a kernel set (and fill it in with provided kernels)
   kernel_set<double> kernels({"sum", "div", "mul","diff"});
   // We can get a vector with these four basic kernels ...
   auto kernel_vector = kernels();
   // (optional) - ... or add our own custom kernel to the set ...
   basic_set.push_back(f);
   // (optional) - ... and get a vector with all five kernels
   auto kernel_vector2 = kernels();

---------------------------------------------------------------------------

.. doxygenclass:: dcgp::kernel_set
   :project: dCGP
   :members:
