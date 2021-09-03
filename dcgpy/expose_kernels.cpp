#include <audi/audi.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <memory>
#include <pagmo/threading.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/function.hpp>
#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/s11n.hpp>
#include <dcgp/wrapped_functions_s11n_implement.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"
#include "pybind11_function.hpp"

using namespace dcgp;
using namespace dcgpy;
using namespace audi;
namespace py = pybind11;

namespace dcgp::s11n_names
{

using udf_py_object_double = dcgp::detail::function_inner<py::object, double, const std::vector<double> &>;
using udf_py_object_string = dcgp::detail::function_inner<py::object, std::string, const std::vector<std::string> &>;
using udf_py_object_gdual_d
    = dcgp::detail::function_inner<py::object, audi::gdual_d, const std::vector<audi::gdual_d> &>;
using udf_py_object_gdual_v
    = dcgp::detail::function_inner<py::object, audi::gdual_v, const std::vector<audi::gdual_v> &>;

using udf_py_object_pc_double
    = dcgp::detail::function_inner<py::object, std::vector<double>, const std::vector<double> &,
                                   std::function<std::vector<double>(const std::vector<double> &)>>;

using udf_py_object_pc_gdual_d
    = dcgp::detail::function_inner<py::object, std::vector<gdual_d>, const std::vector<gdual_d> &,
                                   std::function<std::vector<gdual_d>(const std::vector<gdual_d> &)>>;

using udf_py_object_pc_gdual_v
    = dcgp::detail::function_inner<py::object, std::vector<gdual_v>, const std::vector<gdual_v> &,
                                   std::function<std::vector<gdual_v>(const std::vector<gdual_v> &)>>;

} // namespace dcgp::s11n_names

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_double, "udf py::object double")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_double, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_double)

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_gdual_d, "udf py::object gdual_d")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_gdual_d, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_gdual_d)

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_gdual_v, "udf py::object gdual_v")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_gdual_v, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_gdual_v)

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_string, "udf py::object string")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_string, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_string)

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_pc_double, "udf py::object pc double")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_pc_double, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_pc_double)

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_pc_gdual_d, "udf py::object pc gdual_d")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_pc_gdual_d, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_pc_gdual_d)

BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_py_object_pc_gdual_v, "udf py::object pc gdual_v")
BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_py_object_pc_gdual_v, boost::serialization::track_never)
BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_py_object_pc_gdual_v)

namespace dcgpy
{

// Serialization helpers for kernel.
template <typename T>
py::tuple kernel_pickle_getstate(const kernel<T> &k)
{
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oarchive(oss);
        oarchive << k;
    }
    auto s = oss.str();
    return py::make_tuple(py::bytes(s.data(), boost::numeric_cast<py::size_t>(s.size())));
}

template <typename T>
kernel<T> kernel_pickle_setstate(py::tuple state)
{
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, ("the state tuple passed for kernel deserialization "
                                    "must have 1 element, but instead it has "
                                    + std::to_string(py::len(state)) + " element(s)")
                                       .c_str());
    }

    auto ptr = PyBytes_AsString(state[0].ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "a bytes object is needed to deserialize a kernel");
    }

    std::istringstream iss;
    iss.str(std::string(ptr, ptr + py::len(state[0])));
    kernel<T> k;
    {
        boost::archive::binary_iarchive iarchive(iss);
        iarchive >> k;
    }

    return k;
}

// Serialization helpers for kernel_set.
template <typename T>
py::tuple kernel_set_pickle_getstate(const kernel_set<T> &ks)
{
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oarchive(oss);
        oarchive << ks;
    }
    auto s = oss.str();
    return py::make_tuple(py::bytes(s.data(), boost::numeric_cast<py::size_t>(s.size())));
}

template <typename T>
kernel_set<T> kernel_set_pickle_setstate(py::tuple state)
{
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, ("the state tuple passed for the kernel set deserialization "
                                    "must have 1 element, but instead it has "
                                    + std::to_string(py::len(state)) + " element(s)")
                                       .c_str());
    }

    auto ptr = PyBytes_AsString(state[0].ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "a bytes object is needed to deserialize a kernel set");
    }

    std::istringstream iss;
    iss.str(std::string(ptr, ptr + py::len(state[0])));
    kernel_set<T> ks;
    {
        boost::archive::binary_iarchive iarchive(iss);
        iarchive >> ks;
    }

    return ks;
}

template <typename T>
void expose_kernel(const py::module &m, const std::string &type)
{
    std::string class_name = "kernel_" + type;
    auto ker
        = py::class_<kernel<T>>(m, class_name.c_str(), "The function defining the generic CGP node")
              .def(py::init<>())
              .def(py::init<const py::object &, const py::object &, const std::string &>(), py::arg("callable_f"),
                   py::arg("callable_s"), py::arg("name"), kernel_init_doc(type).c_str())
              .def("__repr__",
                   [](const kernel<T> &instance) -> std::string {
                       std::ostringstream oss;
                       oss << instance;
                       return oss.str();
                   })
              .def(
                  "__call__", [](const kernel<T> &instance, const std::vector<T> &v) { return instance(v); },
                  "Call operator from values")
              .def(
                  "__call__", [](const kernel<T> &instance, const std::vector<std::string> &v) { return instance(v); },
                  "Call operator from strings")

              .def(py::pickle(&dcgpy::kernel_pickle_getstate<T>, &dcgpy::kernel_pickle_setstate<T>));
}

template <typename T>
kernel<T> wrap_operator(const kernel_set<T> &ks, typename std::vector<dcgp::kernel<T>>::size_type idx)
{
    return ks[idx];
}

template <typename T>
void expose_kernel_set(const py::module &m, std::string type)
{
    std::string class_name = "kernel_set_" + type;
    auto ker_set = py::class_<kernel_set<T>>(m, class_name.c_str(),
                                             "Helper to construct a set of kernel functions from their common name");
    ker_set.def(py::init<>())
        .def(py::init<std::vector<std::string>>(), py::arg("kernels"), kernel_set_init_doc(type).c_str())
        .def("__call__", &kernel_set<T>::operator())
        .def("__repr__",
             [](const kernel_set<T> &ks) -> std::string {
                 std::ostringstream oss;
                 oss << ks;
                 return oss.str();
             })
        .def(
            "push_back", [](kernel_set<T> &ks, std::string v) { ks.push_back(v); },
            kernel_set_push_back_str_doc().c_str(), py::arg("kernel_name"))
        .def(
            "push_back", [](kernel_set<T> &ks, const kernel<T> &v) { ks.push_back(v); },
            kernel_set_push_back_ker_doc(type).c_str(), py::arg("kernel"))
        .def("__getitem__", &wrap_operator<T>)
        .def(py::pickle(&dcgpy::kernel_set_pickle_getstate<T>, &dcgpy::kernel_set_pickle_setstate<T>));
}
//
void expose_kernels(const py::module &m)
{
    // double
    expose_kernel<double>(m, "double");
    expose_kernel_set<double>(m, "double");
    //// gdual_d
    expose_kernel<gdual_d>(m, "gdual_double");
    expose_kernel_set<gdual_d>(m, "gdual_double");
    //
    //// gdual_v
    expose_kernel<gdual_v>(m, "gdual_vdouble");
    expose_kernel_set<gdual_v>(m, "gdual_vdouble");
};
} // namespace dcgpy