#ifndef PYDCGP_COMMON_UTILS_HPP
#define PYDCGP_COMMON_UTILS_HPP

#include <boost/numeric/conversion/cast.hpp>
#include <pybind11/pybind11.h>

#include "python_includes.hpp"

namespace py = pybind11;

// Throw a Python exception of type "type" with associated
// error message "msg".
void py_throw(PyObject *type, const char *msg)
{
    PyErr_SetString(type, msg);
    throw py::error_already_set();
}

namespace dcgpy
{

// Perform a deep copy of input object o.
py::object deepcopy(const py::object &o)
{
    return py::module::import("copy").attr("deepcopy")(o);
}

inline std::vector<char> object_to_vchar(const py::object &o)
{
    // This will dump to a bytes object.
    py::object tmp = py::module::import("cloudpickle").attr("dumps")(o);
    // This gives a null-terminated char * to the internal
    // content of the bytes object.
    auto ptr = PyBytes_AsString(tmp.ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "the serialization backend's dumps() function did not return a bytes object");
    }
    // NOTE: this will be the length of the bytes object *without* the terminator.
    const auto size = len(tmp);
    // NOTE: we store as char here because that's what is returned by the CPython function.
    // From Python it seems like these are unsigned chars, but this should not concern us.
    return std::vector<char>(ptr, ptr + size);
}

inline py::object vchar_to_object(const std::vector<char> &v)
{
    auto b = py::bytes(v.data(), boost::numeric_cast<Py_ssize_t>(v.size()));
    return py::module::import("cloudpickle").attr("loads")(b);
}

} // namespace dcgpy

#endif
