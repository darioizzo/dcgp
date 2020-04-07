#ifndef DCGPY_EXPOSE_SYMBOLIC_REGRESSION_HPP
#define DCGPY_EXPOSE_SYMBOLIC_REGRESSION_HPP

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace dcgpy
{
    void expose_symbolic_regression(py::module &);
}

#endif