#include<vector>
#include<string>
#include <functional> //std::function

#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"
#include "pybind11/include/pybind11/functional.h"
#include "../include/basis_function.hpp"

namespace py = pybind11;
using namespace dcgp;

using fun_type = std::function<double(const std::vector<double>&)>;
using fun_print_type = std::function<std::string(const std::vector<std::string>&)>;

PYBIND11_PLUGIN(_core) {
    py::module m("_core", "dCGPy's core module");

    py::class_<basis_function<double>>(m, "kernel")
    .def("__init__",
        [](basis_function<double> &instance, const py::object &obj1, const py::object &obj2, const std::string &name)
            {
                fun_type my_function = [obj1](const std::vector<double>& x) { return obj1(x);};
                fun_print_type my_print_function = [obj2](const std::vector<std::string>& x) { return obj2(x);};
                new (&instance) basis_function<double>(my_function, my_print_function, name);
            }
        )
    .def("__call__",
        [](basis_function<double> &instance, const std::vector<double> in)
        {
            return instance(in);
        }
    );

    return m.ptr();
}
