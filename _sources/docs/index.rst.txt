.. docs tab

Documentation
=============

The documentation is split into two parts. The C++ documentation and the Python documentation. 
We try to keep the APIs as similar as possible, but since templates are not available in python, 
we expose an arbitrary number of template instantiations indicating in the class name
after an underscore, the type. For example, the Python class *expression_double* corresponds, 
to the C++ class *expression<double>*.

-------------------------------------------------

.. toctree::
  :maxdepth: 1

  cpp/index
  python/index
