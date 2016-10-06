#ifndef PYAUDI_COMMON_UTILS_HPP
#define PYAUDI_COMMON_UTILS_HPP

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>



// A throwing macro similar to pagmo_throw, only for Python. This will set the global
// error string of Python to "msg", the exception type to "type", and then invoke the Boost
// Python function to raise the Python exception.
#define dcgpy_throw(type,msg) \
PyErr_SetString(type,msg); \
boost::python::throw_error_already_set(); \
throw

namespace bp = boost::python;

namespace dcgpy{
// Converts a C++ vector to a python list
template <typename T>
inline bp::list v_to_l(std::vector<T> vector) {
    bp::list list;
    for (auto iter = vector.begin(); iter != vector.end(); ++iter) {
        list.append(*iter);
    }
    return list;
}

// Converts a python iterable to an std::vector
template<typename T>
inline std::vector<T> l_to_v(const bp::object& iterable)
{
    bp::stl_input_iterator<T> begin(iterable), end;
    return std::vector<T>(begin, end);
}


}

#endif
