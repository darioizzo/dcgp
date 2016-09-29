#include<vector>
#include<string>
#include<sstream>
#include <functional> //std::function

#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"
#include "pybind11/include/pybind11/functional.h"
#include "../include/basis_function.hpp"
#include "../include/function_set.hpp"
#include "../include/expression.hpp"

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
                fun_type my_function = [obj1](const std::vector<double>& x) { return obj1(x).cast<double>();};
                fun_print_type my_print_function = [obj2](const std::vector<std::string>& x) { return obj2(x).cast<std::string>();};
                new (&instance) basis_function<double>(my_function, my_print_function, name);
            }
        )
    .def("__call__",
        [](basis_function<double> &instance, const std::vector<double> in)
        {
            return instance(in);
        }
    )
    .def("__call__",
        [](basis_function<double> &instance, const std::vector<std::string> in)
        {
            return instance(in);
        }
    )
    .def("__repr__",
        [](const basis_function<double> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    );

    py::class_<function_set<double>>(m, "function_set")
    .def(py::init<>())
    .def(py::init<const std::vector<std::string>&>())
    .def("__call__",
        [](function_set<double> &instance)
        {
            return instance();
        }
    );

    py::class_<expression<double>>(m, "expression")
    .def(py::init<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, std::vector<basis_function<double>>, unsigned int>())
    .def("__repr__",
        [](const expression<double> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    )
    .def("__call__",
        [](const expression<double> &instance, std::vector<double> in)
        {
            return instance(in);
        }
    )
    .def("__call__",
        [](const expression<double> &instance, std::vector<std::string> in)
        {
            return instance(in);
        }
    )
    .def("set", &expression<double>::set, "Sets the expression chromosme")
    .def("get", &expression<double>::get, "Gets the expression chromosme")
    .def("get_lb_", &expression<double>::get_lb, "Gets the lower bounds of the chromosome")
    .def("get_ub", &expression<double>::get_ub, "Gets the upper bounds of the chromosome")
    .def("get_active_genes", &expression<double>::get_active_genes, "Gets the idx of the active genes in the current chromosome (numbering is from 0)")
    .def("get_active_nodes", &expression<double>::get_active_nodes, "Gets the idx of the active nodes in the current chromosome.")
    .def("get_n", &expression<double>::get_n, "Gets the number of inputs of the c_CGP expression")
    .def("get_m", &expression<double>::get_ub, "Gets the number of outputs of the c_CGP expression")
    .def("get_f", &expression<double>::get_f, "Gets the kernel functions")
    .def("mutate", (void (expression<double>::*)(unsigned int)) &expression<double>::mutate, "Mutates the selected gene within its allowed bounds.")
    .def("mutate", (void (expression<double>::*)(std::vector<unsigned int>)) &expression<double>::mutate, "Mutates the selected genes within its allowed bounds.")
    .def("mutate_active", &expression<double>::mutate_active, "Mutates N randomly selected active genes within its allowed bounds.")
    .def("mutate_active_cgene", &expression<double>::mutate_active_cgene, "Mutates exactly one randomly selected active connection genes within its allowed bounds.")
    .def("mutate_ogene", &expression<double>::mutate_ogene, "Mutates exactly one of the output genes within its allowed bounds.")
    .def("mutate_active_fgene", &expression<double>::mutate_active_fgene, "Mutates exactly one randomly selected active function genes within its allowed bounds.");

    return m.ptr();
}
