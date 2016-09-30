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
#include "docstrings.hpp"

namespace py = pybind11;
using namespace dcgp;
using namespace dcgpy;

template <typename T>
void expose_basis_function(const py::module &m, std::string type)
{
    std::string class_name = "kernel_" + type;
    py::class_<basis_function<T>>(m, class_name.c_str())
    .def("__init__",
        [](basis_function<T> &instance, const py::object &obj1, const py::object &obj2, const std::string &name)
            {
                std::function<T(const std::vector<T>&)> my_function = [obj1](const  std::vector<T>& x) { return obj1(x).template cast<T>();};
                std::function<std::string(const std::vector<std::string>&)> my_print_function = [obj2](const std::vector<std::string>& x) { return obj2(x).cast<std::string>();};
                new (&instance) basis_function<T>(my_function, my_print_function, name);
            },
        basis_function_init_doc(type).c_str(),
        py::arg("callable_f"),
        py::arg("callable_s"),
        py::arg("name")
        )
    .def("__call__",
        [](basis_function<T> &instance, const std::vector<T> &in)
        {
            return instance(in);
        }
    )
    .def("__call__",
        [](basis_function<T> &instance, const std::vector<std::string> &in)
        {
            return instance(in);
        }
    )
    .def("__repr__",
        [](const basis_function<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    );
}

template <typename T>
void expose_function_set(const py::module &m, std::string type)
{
    std::string class_name = "function_set_" + type;
    py::class_<function_set<T>>(m, class_name.c_str())
    .def(py::init<>())
    .def(py::init<const std::vector<std::string>&>(),
        function_set_init_doc(type).c_str(),
        py::arg("kernels")
    )
    .def("__call__",
        [](function_set<T> &instance)
        {
            return instance();
        }
    );
}

template <typename T>
void expose_expression(const py::module &m, std::string type)
{
    std::string class_name = "expression_" + type;
    py::class_<expression<T>>(m, class_name.c_str())
    .def(py::init<unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, std::vector<basis_function<T>>, unsigned int>(),
        expression_init_doc(type).c_str(),
        py::arg("in"),
        py::arg("out"),
        py::arg("rows"),
        py::arg("columns"),
        py::arg("levels-back"),
        py::arg("arity"),
        py::arg("kernels"),
        py::arg("seed")
    )
    .def("__repr__",
        [](const expression<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    )
    .def("__call__",
        [](const expression<T> &instance, const std::vector<T> &in)
        {
            return instance(in);
        }
    )
    .def("__call__",
        [](const expression<T> &instance, const std::vector<std::string> &in)
        {
            return instance(in);
        }
    )
    .def("set", &expression<T>::set, "Sets the expression chromosme", py::arg("chromosme"))
    .def("get", &expression<T>::get, "Gets the expression chromosme")
    .def("get_lb_", &expression<T>::get_lb, "Gets the lower bounds of the chromosome")
    .def("get_ub", &expression<T>::get_ub, "Gets the upper bounds of the chromosome")
    .def("get_active_genes", &expression<T>::get_active_genes, "Gets the idx of the active genes in the current chromosome (numbering is from 0)")
    .def("get_active_nodes", &expression<T>::get_active_nodes, "Gets the idx of the active nodes in the current chromosome.")
    .def("get_n", &expression<T>::get_n, "Gets the number of inputs of the c_CGP expression")
    .def("get_m", &expression<T>::get_ub, "Gets the number of outputs of the c_CGP expression")
    .def("get_f", &expression<T>::get_f, "Gets the kernel functions")
    .def("mutate", (void (expression<T>::*)(unsigned int)) &expression<T>::mutate, "Mutates the selected gene within its allowed bounds.", py::arg("idx"))
    .def("mutate", (void (expression<T>::*)(std::vector<unsigned int>)) &expression<T>::mutate, "Mutates the selected genes within its allowed bounds.", py::arg("idxs"))
    .def("mutate_active", &expression<T>::mutate_active, "Mutates N randomly selected active genes within its allowed bounds.")
    .def("mutate_active_cgene", &expression<T>::mutate_active_cgene, "Mutates exactly one randomly selected active connection genes within its allowed bounds.")
    .def("mutate_ogene", &expression<T>::mutate_ogene, "Mutates exactly one of the output genes within its allowed bounds.")
    .def("mutate_active_fgene", &expression<T>::mutate_active_fgene, "Mutates exactly one randomly selected active function genes within its allowed bounds.");

}

PYBIND11_PLUGIN(_core) {
    py::module m("dcgpy", "d-cgpy's core module");

    expose_basis_function<double>(m, "double");
    expose_function_set<double>(m, "double");
    expose_expression<double>(m, "double");
    expose_basis_function<gdual_d>(m, "gdual_d");
    expose_function_set<gdual_d>(m, "gdual_d");
    expose_expression<gdual_d>(m, "gdual_d");

    return m.ptr();
}
