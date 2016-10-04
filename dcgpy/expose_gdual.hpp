#ifndef DCGP_EXPOSE_GDUAL_H
#define DCGP_EXPOSE_GDUAL_H

#include <boost/lexical_cast.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <string>
#include <sstream> //ostringstream, stringstream
#include <stdexcept> // stringstream
#include <vector>

#include <audi/audi.hpp>

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wpedantic"
    #pragma GCC diagnostic ignored "-Wshadow"
    #pragma GCC diagnostic ignored "-Wsign-conversion"
    #pragma GCC diagnostic ignored "-Wdeprecated"
#endif

#include "pybind11/include/pybind11/operators.h"
#include "pybind11/include/pybind11/pybind11.h"
#include "pybind11/include/pybind11/stl.h"
#include "pybind11/include/pybind11/cast.h"

#if defined(__clang__) || defined(__GNUC__)
    #pragma GCC diagnostic pop
#endif

using namespace audi;
namespace py = pybind11;

namespace pyaudi {

template<typename T>
auto expose_gdual(py::module &m, std::string type)
{
     return py::class_<gdual<T>>(m,("gdual_"+type).c_str())
    .def(py::init<>())
    .def(py::init<const gdual<T> &>())
    .def(py::init<T>())
    .def(py::init<T, const std::string &, unsigned int>())
    .def("__repr__",[](const gdual<T> &g) -> std::string {
        std::ostringstream oss;
        oss << g;
        return oss.str();
    })
    .def("_repr_latex_",[](const gdual<T> &g) -> std::string {
        std::ostringstream oss;
        g._poly().print_tex(oss);
        auto retval = oss.str();
        retval += std::string("+\\mathcal{O}\\left(")
            + boost::lexical_cast<std::string>(g.get_order() + 1) +  "\\right) \\]";
        return std::string("\\[ ") + retval;
    })
    .def("__getstate__", [](const gdual<T> &p) {
        // Returns a tuple that contains the string
        // representation of a gdual<T> as obtained
        // from the boost serialization library
        std::stringstream ss;
        boost::archive::text_oarchive oa(ss);
        oa << p;
        return py::make_tuple(ss.str());
    })
    .def("__setstate__", [](gdual<T> &p, py::tuple t) {
        if (t.size() != 1)
            throw std::runtime_error("Invalid state!");
        // Invoke the default constructor.
        new (&p) gdual<T>;
        // Reconstruct the gdual<T>
        std::stringstream ss(t[0].cast<std::string>());
        boost::archive::text_iarchive ia(ss);
        ia >> p;
    })
    .def_property_readonly("symbol_set",&gdual<T>::get_symbol_set, "The list of symbols in the polynomial")
    .def_property_readonly("symbol_set_size",&gdual<T>::get_symbol_set_size)
    .def_property_readonly("degree",&gdual<T>::degree, "polynomial degree (<= order)")
    .def_property_readonly("order",&gdual<T>::get_order, "truncation order (>= degree)")
    .def_property_readonly("constant_cf",&gdual<T>::constant_cf, "Constant term of the polynomial")
    .def("extend_symbol_set", &gdual<T>::extend_symbol_set, "Extends the symbol set")
    .def("integrate", &gdual<T>::template integrate<>, "Integrate with respect to argument")
    .def("partial", &gdual<T>::template partial<>, "Partial derivative with respect to argument")
    .def("evaluate",[](const gdual<T> &g, const std::map< std::string, double> &dict) {return g.evaluate(std::unordered_map< std::string, double>(dict.begin(), dict.end()));} , "Evaluates the Taylor polynomial")
    .def("find_cf", [](const gdual<T> &g, const std::vector<int> &v) {
        return g.find_cf(v);
    },"Find the coefficient of the Taylor expansion")
    .def("get_derivative", [](const gdual<T> &g, const std::vector<int> &v) {
        return g.get_derivative(v);
    },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
    .def("get_derivative", [](const gdual<T> &g, const std::unordered_map<std::string, unsigned int> &dict) {
        return g.get_derivative(dict);
    },"Finds the derivative (i.e. the coefficient of the Taylor expansion discounted of a factorial factor")
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self * py::self)
    .def(py::self / py::self)
    .def(py::self + double())
    .def(py::self - double())
    .def(py::self * double())
    .def(py::self / double())
    .def(-py::self)
    .def(+py::self)
    .def(double() + py::self)
    .def(double() - py::self)
    .def(double() * py::self)
    .def(double() / py::self)
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def("__pow__",[](const gdual<T> &g, double x) {return pow(g,x);} ,("Exponentiation (gdual_"+type+", double).").c_str())
    .def("__pow__",[](const gdual<T> &base, const gdual<T> &g) {return pow(base,g);} ,("Exponentiation (gdual_"+type+", gdual_"+type+").").c_str())
    .def("__rpow__",[](const gdual<T> &g, double x) {return pow(x,g);} ,("Exponentiation (double, gdual_"+type+").").c_str())
    ;
}
}
#endif
