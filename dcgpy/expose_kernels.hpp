#ifndef DCGPY_EXPOSE_KERNELS_HPP
#define DCGPY_EXPOSE_KERNELS_HPP

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace dcgpy
{
    void expose_kernels(const py::module &);
}

#endif