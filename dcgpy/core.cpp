#include <cmath>
#include <pybind11/pybind11.h>
#include <sstream>

#include <dcgp/problems/symbolic_regression.hpp>

#include <boost/optional.hpp>

#include <tbb/global_control.h>

#include "docstrings.hpp"
#include "expose_expressions.hpp"
#include "expose_kernels.hpp"
#include "expose_symbolic_regression.hpp"


using namespace dcgpy;
namespace py = pybind11;

namespace {
    boost::optional<tbb::global_control> thread_control; 
}

PYBIND11_MODULE(core, m)
{
    py::options options;
    options.disable_function_signatures();

    expose_kernels(m);
    expose_expressions(m);
    expose_symbolic_regression(m);

    // Override the default implementation of the island factory.
    dcgp::details::extract_sr_cpp_py = [](const pagmo::problem &p) -> const dcgp::symbolic_regression * {
        auto py_ptr = p.extract<py::object>();
        if (!py_ptr) {
            return nullptr;
        }
        const dcgp::symbolic_regression *retval;
        try {
            retval = py_ptr->cast<const dcgp::symbolic_regression *>();
        } catch (py::cast_error) {
            retval = nullptr;
        }
        return retval;
    };

    m.def("disable_threading", [](){ thread_control.emplace(tbb::global_control::max_allowed_parallelism, 1); }, disable_threading_doc().c_str());
    m.def("enable_threading", [](){ thread_control.reset(); }, enable_threading_doc().c_str());

} // namespace details
