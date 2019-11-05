// See: https://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
// In every cpp file We need to make sure this is included before everything else,
// with the correct #defines.
#include "python_includes.hpp"
#define PY_ARRAY_UNIQUE_SYMBOL dcgpy_ARRAY_API
#include "numpy.hpp"

#include <audi/audi.hpp>
#include <boost/python.hpp>
#include <pygmo/register_ap.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"
#include "expose_expressions.hpp"
#include "expose_kernels.hpp"
#include "expose_symbolic_regression.hpp"

using namespace dcgpy;
using namespace audi;
namespace bp = boost::python;
namespace pg = pygmo;

BOOST_PYTHON_MODULE(core)
{
    bp::docstring_options doc_options;

    // Init numpy.
    // NOTE: only the second import is strictly necessary. We run a first import from BP
    // because that is the easiest way to detect whether numpy is installed or not (rather
    // than trying to figure out a way to detect it from import_array()).
    try {
        bp::import("numpy.core.multiarray");
    } catch (...) {
        dcgpy::builtin().attr("print")(
            u8"\033[91m====ERROR====\nThe NumPy module could not be imported. "
            u8"Please make sure that NumPy has been correctly installed.\n====ERROR====\033[0m");
        throw std::runtime_error("Import error");
    }
    dcgpy::numpy_import_array();

    doc_options.enable_all();
    doc_options.disable_cpp_signatures();
    doc_options.disable_py_signatures();

    // Registers dcgpy as a pygmo affiliated package, so that upon import will add itself
    // to the cereal serialization dictionary
    pg::register_ap();

    // Expose all expressions
    dcgpy::expose_expressions();
    // Expose all kernels and kernel_sets
    dcgpy::expose_kernels();
    // Expose Symbolic Regression Stuff
    dcgpy::expose_symbolic_regression();

    // Define a cleanup functor to be run when the module is unloaded.
    struct dcgp_cleanup_functor {
        void operator()() const
        {
            std::cout << "Shutting down the thread pool.\n";
            piranha::thread_pool_shutdown<void>();
        }
    };
    // Expose it.
    bp::class_<dcgp_cleanup_functor> cl_c("_dcgp_cleanup_functor", bp::init<>());
    cl_c.def("__call__", &dcgp_cleanup_functor::operator());
    // Register it.
    bp::object atexit_mod = bp::import("atexit");
    atexit_mod.attr("register")(dcgp_cleanup_functor{});
}
