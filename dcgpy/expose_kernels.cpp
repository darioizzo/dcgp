// See: https://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
// In every cpp file We need to make sure this is included before everything else,
// with the correct #defines.
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL dcgpy_ARRAY_API
#include "numpy.hpp"

#include <boost/python.hpp>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"

using namespace dcgp;
using namespace dcgpy;
using namespace audi;
namespace bp = boost::python;

namespace dcgpy
{

template <typename T>
void expose_kernel(const std::string &type)
{
    std::string class_name = "kernel_" + type;
    bp::class_<kernel<T>>(class_name.c_str(), "The function defining the generic CGP node", bp::no_init)
        .def("__init__",
             bp::make_constructor(
                 +[](const bp::object &obj1, const bp::object &obj2, const std::string &name) {
                     std::function<T(const std::vector<T> &)> my_function = [obj1](const std::vector<T> &x) {
                         T in = bp::extract<T>(obj1(v_to_l(x)));
                         return in;
                     };
                     std::function<std::string(const std::vector<std::string> &)> my_print_function
                         = [obj2](const std::vector<std::string> &x) {
                               std::string in = bp::extract<std::string>(obj2(v_to_l(x)));
                               return in;
                           };
                     return ::new kernel<T>(my_function, my_print_function, name);
                 },
                 bp::default_call_policies(), (bp::arg("callable_f"), bp::arg("callable_s"), bp::arg("name"))),
             kernel_init_doc(type).c_str())
        .def(
            "__call__",
            +[](kernel<T> &instance, const bp::object &in) {
                try {
                    auto v = l_to_v<T>(in);
                    return bp::object(instance(v));
                } catch (...) {
                    PyErr_Clear();
                    auto v = l_to_v<std::string>(in);
                    return bp::object(instance(v));
                }
            })
        .def(
            "__repr__", +[](const kernel<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            });
    ;
}

template <typename T>
kernel<T> wrap_operator(const kernel_set<T> &ks, typename std::vector<dcgp::kernel<T>>::size_type idx)
{
    return ks[idx];
}

template <typename T>
void expose_kernel_set(std::string type)
{
    std::string class_name = "kernel_set_" + type;
    bp::class_<kernel_set<T>>(class_name.c_str(),
                              "Helper to construct a set of kernel functions from their common name", bp::no_init)
        .def("__init__",
             bp::make_constructor(
                 +[](const bp::object &obj1) {
                     auto a = l_to_v<std::string>(obj1);
                     return ::new kernel_set<T>(a);
                 },
                 bp::default_call_policies(), (bp::arg("kernels"))),
             kernel_set_init_doc(type).c_str())
        .def(
            "__call__", +[](kernel_set<T> &instance) { return v_to_l(instance()); })
        .def(
            "__repr__",
            +[](const kernel_set<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            })
        .def("push_back", (void (kernel_set<T>::*)(std::string)) & kernel_set<T>::push_back,
             kernel_set_push_back_str_doc().c_str(), bp::arg("kernel_name"))
        .def("push_back", (void (kernel_set<T>::*)(const kernel<T> &)) & kernel_set<T>::push_back,
             kernel_set_push_back_ker_doc(type).c_str(), bp::arg("kernel"))
        .def("__getitem__", &wrap_operator<T>);
}

void expose_kernels()
{
    // double
    expose_kernel<double>("double");
    expose_kernel_set<double>("double");

    // gdual_d
    expose_kernel<gdual_d>("gdual_double");
    expose_kernel_set<gdual_d>("gdual_double");

    // gdual_v
    expose_kernel<gdual_v>("gdual_vdouble");
    expose_kernel_set<gdual_v>("gdual_vdouble");
}
} // namespace dcgpy