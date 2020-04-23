#ifndef PYDCGP_COMMON_UTILS_HPP
#define PYDCGP_COMMON_UTILS_HPP

#include <algorithm>
#include <boost/numeric/conversion/cast.hpp>
#include <pagmo/types.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

#include "python_includes.hpp"

namespace py = pybind11;

namespace dcgpy
{
// Throw a Python exception of type "type" with associated
// error message "msg".
[[noreturn]] inline void py_throw(PyObject *type, const char *msg)
{
    PyErr_SetString(type, msg);
    throw py::error_already_set();
}

// Perform a deep copy of input object o.
inline py::object deepcopy(const py::object &o)
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

// Converts a numpy array into a vector
inline std::vector<double> ndarr_to_vector(const py::array_t<double> &a)
{
    // Get a two-dimensional view on the array.
    // If the array is not 2D, this will throw.
    auto r = a.template unchecked<1>();
    auto n = r.shape(0);

    // Prepare the retval.
    std::vector<double> retval(static_cast<std::vector<std::vector<double>>::size_type>(n), 0.);

    // Copy the values
    for (py::ssize_t j = 0u; j < n; ++j) {
        retval[static_cast<std::vector<double>::size_type>(j)] = r(j);
    }

    return retval;
}

// Convert an std::vector<std::vector> into a 2-D numpy array.
inline py::array_t<double> vector_to_ndarr(const std::vector<double> &vv)
{
    // Check data
    auto n = vv.size();
    if (n == 0u) {
        py_throw(PyExc_ValueError, "Empty C++ vectors cannot be converted to Numpy arrays.");
    }

    // Create the output array, of shape n
    py::array_t<double> retval(n);

    // Get a mutable view into it and copy the data from vv.
    auto r = retval.mutable_unchecked<1>();
    for (decltype(n) i = 0u; i < n; ++i) {
        r(static_cast<py::ssize_t>(i)) = vv[i];
    }
    return retval;
}

// Converts a numpy array into a vector of vectors.
inline std::vector<std::vector<double>> ndarr_to_vvector(const py::array_t<double> &a)
{
    // Get a two-dimensional view on the array.
    // If the array is not 2D, this will throw.
    auto r = a.template unchecked<2>();
    auto n = r.shape(0);
    auto m = r.shape(1);

    // Prepare the retval.
    std::vector<std::vector<double>> retval(static_cast<std::vector<std::vector<double>>::size_type>(n),
                                            std::vector<double>(static_cast<std::vector<double>::size_type>(m), 0.));

    // Copy the values
    for (py::ssize_t i = 0u; i < n; ++i) {
        for (py::ssize_t j = 0u; j < m; ++j) {
            retval[static_cast<std::vector<std::vector<double>>::size_type>(i)]
                  [static_cast<std::vector<double>::size_type>(j)]
                = r(i, j);
        }
    }
    return retval;
}

// Convert an std::vector<std::vector> into a 2-D numpy array.
inline py::array_t<double> vvector_to_ndarr(const std::vector<std::vector<double>> &vv)
{
    // Check data
    auto n = vv.size();
    if (n == 0u) {
        py_throw(PyExc_ValueError, "Empty C++ vectors cannot be converted to Numpy arrays.");
    }
    auto m = vv[0].size();
    if (m == 0u) {
        py_throw(PyExc_ValueError, "Empty C++ vectors cannot be converted to Numpy arrays.");
    }
    if (std::any_of(vv.cbegin(), vv.cend(), [m](auto v) { return v.size() != m; })) {
        py_throw(PyExc_ValueError, "Data is malformed for a Numpy array conversion (i.e. sizes are not consistent)");
    }

    // Create the output array, of shape n x m.
    py::array_t<double> retval({boost::numeric_cast<py::ssize_t>(n), boost::numeric_cast<py::ssize_t>(m)});

    // Get a mutable view into it and copy the data from vv.
    auto r = retval.mutable_unchecked<2>();
    for (decltype(n) i = 0u; i < n; ++i) {
        for (decltype(m) j = 0u; j < m; ++j) {
            r(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j)) = vv[i][j];
        }
    }

    return retval;
}

// Convert a sparsity pattern into a numpy array.
inline py::array_t<pagmo::vector_double::size_type> sp_to_ndarr(const pagmo::sparsity_pattern &sp)
{
    // Create the output array, of shape n x 2.
    py::array_t<pagmo::vector_double::size_type> retval({boost::numeric_cast<py::ssize_t>(sp.size()), py::ssize_t(2)});

    // Get a mutable view into it and copy the data from sp.
    auto r = retval.mutable_unchecked<2>();
    for (decltype(sp.size()) i = 0; i < sp.size(); ++i) {
        r(static_cast<py::ssize_t>(i), 0) = sp[i].first;
        r(static_cast<py::ssize_t>(i), 1) = sp[i].second;
    }

    return retval;
}

// Utils to expose algo log.
template <typename Algo>
inline py::list generic_log_getter(const Algo &a)
{
    py::list retval;
    for (const auto &t : a.get_log()) {
        retval.append(t);
    }
    return retval;
}

template <typename Algo>
inline void expose_algo_log(py::class_<Algo> &algo_class, const char *doc)
{
    algo_class.def("get_log", &generic_log_getter<Algo>, doc);
}

// Serialization support
template <typename UDX>
inline py::tuple udx_pickle_getstate(const UDX &udx)
{
    // The idea here is that first we extract a char array
    // into which a has been serialized, then we turn
    // this object into a Python bytes object and return that.
    std::ostringstream oss;
    {
        boost::archive::binary_oarchive oarchive(oss);
        oarchive << udx;
    }
    auto s = oss.str();
    return py::make_tuple(py::bytes(s.data(), boost::numeric_cast<py::size_t>(s.size())));
}

template <typename UDX>
inline UDX udx_pickle_setstate(py::tuple state)
{
    // Similarly, first we extract a bytes object from the Python state,
    // and then we build a C++ string from it. The string is then used
    // to deserialized the object.
    if (py::len(state) != 1) {
        py_throw(PyExc_ValueError, ("the state tuple passed for udp/uda deserialization "
                                    "must have 1 element, but instead it has "
                                    + std::to_string(py::len(state)) + " element(s)")
                                       .c_str());
    }

    auto ptr = PyBytes_AsString(state[0].ptr());
    if (!ptr) {
        py_throw(PyExc_TypeError, "a bytes object is needed to deserialize a problem / algorithm");
    }

    std::istringstream iss;
    iss.str(std::string(ptr, ptr + py::len(state[0])));
    UDX udx;
    {
        boost::archive::binary_iarchive iarchive(iss);
        iarchive >> udx;
    }

    return udx;
}

} // namespace dcgpy

#endif
