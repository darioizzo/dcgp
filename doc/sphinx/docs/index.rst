.. docs tab

Documentation
=============

The documentation is split into two parts. The C++ documentation and the Python documentation. We try to keep the APIs as similar as possible, but since temmplates are not available in python, we expose
an arbitrary number of template instantiations indoicasting after an underscore the type. For example, the class *expression_double* corresponds, in Python, to the C++ class *expression<double>*.

-------------------------------------------------

.. toctree::
  :maxdepth: 3

  cpp_docs

------------------------------------------------

.. toctree::
  :maxdepth: 2

  python_docs
