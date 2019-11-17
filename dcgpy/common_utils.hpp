#ifndef PYDCGP_COMMON_UTILS_HPP
#define PYDCGP_COMMON_UTILS_HPP

#include <boost/numeric/conversion/cast.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

#include <string>
#include <vector>

#include "numpy.hpp"

// A throwing macro similar to pagmo_throw, only for Python. This will set the global
// error string of Python to "msg", the exception type to "type", and then invoke the Boost
// Python function to raise the Python exception.
#define dcgpy_throw(type, msg)                                                                                         \
    PyErr_SetString(type, msg);                                                                                        \
    boost::python::throw_error_already_set();                                                                          \
    throw

namespace bp = boost::python;

// Map C++ types to NPY_ types.
template <typename T>
struct cpp_npy {
};

#define DCGPY_CPP_NPY(from, to)                                                                                        \
    template <>                                                                                                        \
    struct cpp_npy<from> {                                                                                             \
        static constexpr auto value = to;                                                                              \
    };

// We only need the types below at the moment.
DCGPY_CPP_NPY(unsigned char, NPY_UBYTE)
DCGPY_CPP_NPY(unsigned short, NPY_USHORT)
DCGPY_CPP_NPY(unsigned, NPY_UINT)
DCGPY_CPP_NPY(unsigned long, NPY_ULONG)
DCGPY_CPP_NPY(unsigned long long, NPY_ULONGLONG)
DCGPY_CPP_NPY(signed char, NPY_BYTE)
DCGPY_CPP_NPY(short, NPY_SHORT)
DCGPY_CPP_NPY(int, NPY_INT)
DCGPY_CPP_NPY(long, NPY_LONG)
DCGPY_CPP_NPY(long long, NPY_LONGLONG)
DCGPY_CPP_NPY(float, NPY_FLOAT)
DCGPY_CPP_NPY(double, NPY_DOUBLE)

#undef DCGPY_CPP_NPY

namespace dcgpy
{

// Import and return the builtin module.
inline bp::object builtin()
{
#if PY_MAJOR_VERSION < 3
    return bp::import("__builtin__");
#else
    return bp::import("builtins");
#endif
}

// isinstance wrapper.
inline bool isinstance(const bp::object &o, const bp::object &t)
{
    return bp::extract<bool>(builtin().attr("isinstance")(o, t));
}

// String representation of an object.
inline std::string str(const bp::object &o)
{
    return bp::extract<std::string>(builtin().attr("str")(o));
}

// Get the type of an object.
inline bp::object type(const bp::object &o)
{
    return builtin().attr("type")(o);
}

// Converts a C++ vector to a python list
template <typename T>
inline bp::list v_to_l(std::vector<T> vector)
{
    bp::list list;
    for (auto iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}

// Converts a python iterable to an std::vector
template <typename T>
inline std::vector<T> l_to_v(const bp::object &iterable)
{
    bp::stl_input_iterator<T> begin(iterable), end;
    return std::vector<T>(begin, end);
}

inline std::vector<double> ad_to_vd(PyArrayObject *o)
{
    assert(PyArray_TYPE(o) == NPY_DOUBLE);
    using size_type = std::vector<double>::size_type;
    if (!PyArray_ISCARRAY_RO(o)) {
        throw std::runtime_error("cannot convert NumPy array to a vector of doubles: "
                                 "data must be C-style contiguous, aligned, and in machine byte-order");
    }
    if (PyArray_NDIM(o) != 1) {
        throw std::runtime_error("cannot convert NumPy array to a vector of doubles: "
                                 "the array must be unidimensional, but the dimension is "
                                 + std::to_string(PyArray_NDIM(o)) + " instead");
    }
    if (PyArray_STRIDES(o)[0] != sizeof(double)) {
        throw std::runtime_error("cannot convert NumPy array to a vector of doubles: "
                                 "the stride value must be "
                                 + std::to_string(sizeof(double)));
    }
    if (PyArray_ITEMSIZE(o) != sizeof(double)) {
        throw std::runtime_error("cannot convert NumPy array to a vector of doubles: "
                                 "the size of the scalar type must be "
                                 + std::to_string(sizeof(double)));
    }
    // NOTE: not sure if this special casing is needed. We make sure
    // the array contains something in order to avoid messing around
    // with a potentially null pointer in the array.
    const auto size = boost::numeric_cast<size_type>(PyArray_SHAPE(o)[0]);
    if (size) {
        auto data = static_cast<double *>(PyArray_DATA(o));
        return std::vector<double>(data, data + size);
    }
    return std::vector<double>{};
}

// Converts an arbitrary python iterable to a vector<T>. Works for NumpyArray too, but not efficient.
template <typename T>
inline std::vector<T> to_v(const bp::object &o)
{
    bp::stl_input_iterator<T> begin(o), end;
    return std::vector<T>(begin, end);
}

// Template specialization for the case of doubles where we distinguish the case of NumpyArray
template <>
inline std::vector<double> to_v<double>(const bp::object &o)
{
    bp::object a = bp::import("numpy").attr("ndarray");
    if (isinstance(o, a)) {
        // NOTE: the idea here is that we want to be able to convert
        // from a NumPy array of types other than double. This is useful
        // because one can then create arrays of ints and have them converted
        // on the fly (e.g., for the bounds). If the array is already a
        // double-precision array, this function should not do any copy.
        auto n = PyArray_FROM_OTF(o.ptr(), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (!n) {
            bp::throw_error_already_set();
        }
        return ad_to_vd(reinterpret_cast<PyArrayObject *>(bp::object(bp::handle<>(n)).ptr()));
    }
    // If o is not a numpy array, just try to iterate over it and extract doubles.
    bp::stl_input_iterator<double> begin(o), end;
    return std::vector<double>(begin, end);
}

// Converts a 2-D Numpy array of doubles to a vector<vector<double>>.
inline std::vector<std::vector<double>> ad_to_vvd(PyArrayObject *o)
{
    using size_type = std::vector<std::vector<double>>::size_type;
    if (!PyArray_ISCARRAY_RO(o)) {
        throw std::runtime_error(
            "cannot convert NumPy array to a vector of vector_double: data must be C-style contiguous, "
            "aligned, and in machine byte-order");
    }
    if (PyArray_NDIM(o) != 2) {
        throw std::invalid_argument(
            "cannot convert NumPy array to a vector of vector_double: the array must be 2-dimensional");
    }
    if (PyArray_TYPE(o) != NPY_DOUBLE) {
        throw std::invalid_argument(
            "cannot convert NumPy array to a vector of vector_double: the scalar type must be double");
    }
    if (PyArray_ITEMSIZE(o) != sizeof(double)) {
        throw std::runtime_error(
            "cannot convert NumPy array to a vector of vector_double:  the size of the scalar type must be "
            + std::to_string(sizeof(double)));
    }
    const auto size = boost::numeric_cast<size_type>(PyArray_SHAPE(o)[0]);
    std::vector<std::vector<double>> retval;
    if (size) {
        auto data = static_cast<double *>(PyArray_DATA(o));
        const auto ssize = PyArray_SHAPE(o)[1];
        for (size_type i = 0u; i < size; ++i, data += ssize) {
            retval.push_back(std::vector<double>(data, data + ssize));
        }
    }
    return retval;
}

// Convert an arbitrary Python object to a vector of vectors.
template <typename T>
inline std::vector<std::vector<T>> to_vv(const bp::object &o)
{
    bp::stl_input_iterator<bp::object> begin(o), end;
    std::vector<std::vector<T>> retval;
    for (; begin != end; ++begin) {
        retval.push_back(to_v<T>(*begin));
    }
    return retval;
}

// Template specialization for the case of doubles where we distinguish the case of NumpyArray.
template <>
inline std::vector<std::vector<double>> to_vv<double>(const bp::object &o)
{
    bp::object l = builtin().attr("list");
    bp::object a = bp::import("numpy").attr("ndarray");
    if (isinstance(o, l)) {
        bp::stl_input_iterator<bp::object> begin(o), end;
        std::vector<std::vector<double>> retval;
        for (; begin != end; ++begin) {
            retval.push_back(to_v<double>(*begin));
        }
        return retval;
    } else if (isinstance(o, a)) {
        auto n = PyArray_FROM_OTF(o.ptr(), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
        if (!n) {
            bp::throw_error_already_set();
        }
        return ad_to_vvd(reinterpret_cast<PyArrayObject *>(bp::object(bp::handle<>(n)).ptr()));
    }
    throw std::invalid_argument(
        "cannot convert the type '" + str(type(o))
        + "' to a vector of vector_double: only lists of doubles and NumPy arrays of doubles are supported");
}

// Convert a vector of arithmetic types into a 1D numpy array.
template <typename T>
inline bp::object vector_to_ndarr(const std::vector<T> &v)
{
    // The dimensions of the array to be created.
    npy_intp dims[] = {boost::numeric_cast<npy_intp>(v.size())};
    // Attempt creating the array.
    PyObject *ret = PyArray_SimpleNew(1, dims, cpp_npy<T>::value);
    if (!ret) {
        dcgpy_throw(PyExc_RuntimeError, "couldn't create a NumPy array: the 'PyArray_SimpleNew()' function failed");
    }
    // Hand over to BP for exception-safe behaviour.
    bp::object retval{bp::handle<>(ret)};
    if (v.size()) {
        // Copy over the data.
        std::copy(v.begin(), v.end(), static_cast<T *>(PyArray_DATA(reinterpret_cast<PyArrayObject *>(ret))));
    }
    return retval;
}

// Convert a vector of vectors of arithmetic types into a 2D numpy array.
template <typename T>
inline bp::object vvector_to_ndarr(const std::vector<std::vector<T>> &v)
{
    // The dimensions of the array to be created.
    const auto nrows = v.size();
    const auto ncols = nrows ? v[0].size() : 0u;
    npy_intp dims[] = {boost::numeric_cast<npy_intp>(nrows), boost::numeric_cast<npy_intp>(ncols)};
    // Attempt creating the array.
    PyObject *ret = PyArray_SimpleNew(2, dims, cpp_npy<T>::value);
    if (!ret) {
        dcgpy_throw(PyExc_RuntimeError, "couldn't create a NumPy array: the 'PyArray_SimpleNew()' function failed");
    }
    // Hand over to BP for exception-safe behaviour.
    bp::object retval{bp::handle<>(ret)};
    if (nrows) {
        auto data = static_cast<T *>(PyArray_DATA(reinterpret_cast<PyArrayObject *>(ret)));
        for (const auto &i : v) {
            if (i.size() != ncols) {
                dcgpy_throw(PyExc_ValueError, "cannot convert a vector of vectors to a NumPy 2D array "
                                              "if the vector instances don't have all the same size");
            }
            std::copy(i.begin(), i.end(), data);
            data += ncols;
        }
    }
    return retval;
}

// Wrapper around the CPython function to create a bytes object from raw data.
inline bp::object make_bytes(const char *ptr, Py_ssize_t len)
{
    PyObject *retval;
    if (len) {
        retval = PyBytes_FromStringAndSize(ptr, len);
    } else {
        retval = PyBytes_FromStringAndSize(nullptr, 0);
    }
    if (!retval) {
        dcgpy_throw(PyExc_RuntimeError, "unable to create a bytes object: the 'PyBytes_FromStringAndSize()' "
                                        "function returned NULL");
    }
    return bp::object(bp::handle<>(retval));
}

// Perform a deep copy of input object o.
inline bp::object deepcopy(const bp::object &o)
{
    return bp::import("copy").attr("deepcopy")(o);
}

inline std::vector<char> object_to_vchar(const bp::object &o)
{
    // This will dump to a bytes object.
    bp::object tmp = bp::import("cloudpickle").attr("dumps")(o);
    // This gives a null-terminated char * to the internal
    // content of the bytes object.
    auto ptr = PyBytes_AsString(tmp.ptr());
    if (!ptr) {
        dcgpy_throw(PyExc_TypeError, "the serialization backend's dumps() function did not return a bytes object");
    }
    // NOTE: this will be the length of the bytes object *without* the terminator.
    const auto size = len(tmp);
    // NOTE: we store as char here because that's what is returned by the CPython function.
    // From Python it seems like these are unsigned chars, but this should not concern us.
    return std::vector<char>(ptr, ptr + size);
}

inline bp::object vchar_to_object(const std::vector<char> &v)
{
    auto b = make_bytes(v.data(), boost::numeric_cast<Py_ssize_t>(v.size()));
    return bp::import("cloudpickle").attr("loads")(b);
}

} // namespace dcgpy

#endif
