#ifndef DCGPY_EXPOSE_EXPRESSIONS_HPP
#define DCGPY_EXPOSE_EXPRESSIONS_HPP

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace dcgpy
{
    void expose_expressions(const py::module &);
}

#endif