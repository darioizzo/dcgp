#include <dcgp/dcgp.hpp>
#include <cmath>
#include <pybind11/pybind11.h>

#include <sstream>

//#include "common_utils.hpp"
#include "expose_kernels.hpp"

using namespace dcgpy;
namespace py = pybind11;

PYBIND11_MODULE(core, m)
{
    expose_kernels(m);
}