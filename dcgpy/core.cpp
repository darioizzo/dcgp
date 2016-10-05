#include <audi/audi.hpp>
#include <boost/python.hpp>
#include <vector>
#include <string>
#include <sstream>
#include <functional> //std::function

#include "common_utils.hpp"
#include "../include/basis_function.hpp"
#include "../include/function_set.hpp"
#include "../include/expression.hpp"
#include "docstrings.hpp"

using namespace dcgp;
using namespace dcgpy;
using namespace audi;
namespace bp = boost::python;

template <typename T>
void expose_basis_function(const std::string &type)
{
    std::string class_name = "kernel_" + type;
    bp::class_<basis_function<T>>(class_name.c_str(), bp::no_init)
    .def("__init__", bp::make_constructor(
        +[](const bp::object &obj1, const bp::object &obj2, const std::string &name)
            {
                std::function<T(const std::vector<T>&)> my_function = [obj1](const std::vector<T>& x)
                {
                    T in = bp::extract<T>(obj1(v_to_l(x)));
                    return in;
                };
                std::function<std::string(const std::vector<std::string>&)> my_print_function = [obj2](const std::vector<std::string>& x)
                {
                    std::string in = bp::extract<std::string>(obj2(v_to_l(x)));
                    return in;
                };
                return ::new basis_function<T>(my_function, my_print_function, name);
            },
            bp::default_call_policies(),
            (bp::arg("callable_f"), bp::arg("callable_s"), bp::arg("name"))
            ),
            basis_function_init_doc(type).c_str()
        )
    .def("__call__",
        +[](basis_function<T> &instance, const bp::object &in)
        {
            try {
                auto v = l_to_v<double>(in);
                return bp::object(instance(v));
            } catch (...) {
                auto v = l_to_v<std::string>(in);
                return bp::object(instance(v));
            }
        }
    )
    .def("__repr__",
        +[](const basis_function<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    );
}


BOOST_PYTHON_MODULE(_core)
{
    expose_basis_function<double>("double");
}
