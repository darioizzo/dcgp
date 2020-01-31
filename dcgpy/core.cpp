#include <cmath>
#include <pybind11/pybind11.h>
#include <sstream>

#include <dcgp/dcgp.hpp>

#include "expose_expressions.hpp"
#include "expose_kernels.hpp"

using namespace dcgpy;
namespace py = pybind11;

PYBIND11_MODULE(core, m)
{
    expose_kernels(m);
    expose_expressions(m);
}