// See: https://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
// In every cpp file We need to make sure this is included before everything else,
// with the correct #defines.
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL dcgpy_ARRAY_API
#include "numpy.hpp"

#include <boost/python.hpp>
#include <pygmo/algorithm_exposition_suite.hpp>
#include <pygmo/problem_exposition_suite.hpp>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/gym.hpp>
#include <dcgp/kernel.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"

namespace bp = boost::python;
namespace pg = pygmo;
using namespace dcgp;

using gym_ptr = void (*)(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
template <gym_ptr F>
inline void expose_data_from_the_gym(const std::string &name, const std::string &docstring)
{
    bp::def(
        name.c_str(),
        +[]() {
            std::vector<std::vector<double>> points, labels;
            F(points, labels);
            return bp::make_tuple(dcgpy::vvector_to_ndarr<double>(points), dcgpy::vvector_to_ndarr<double>(labels));
        },
        docstring.c_str());
}

namespace dcgpy
{

void expose_symbolic_regression()
{
    // We expose the UDPs
    pg::expose_problem<dcgp::symbolic_regression>("symbolic_regression", symbolic_regression_doc().c_str())
        .def("__init__",
             bp::make_constructor(
                 +[](const bp::object &points, const bp::object &labels, unsigned rows, unsigned cols,
                     unsigned levels_back, unsigned arity, const bp::object &kernels, unsigned n_eph,
                     bool multi_objective, unsigned parallel_batches) {
                     auto kernels_v = l_to_v<kernel<double>>(kernels);
                     auto vvd_points = to_vv<double>(points);
                     auto vvd_labels = to_vv<double>(labels);
                     return ::new dcgp::symbolic_regression(vvd_points, vvd_labels, rows, cols, levels_back, arity,
                                                            kernels_v, n_eph, multi_objective, parallel_batches);
                 },
                 bp::default_call_policies(),
                 (bp::arg("points"), bp::arg("labels"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"), bp::arg("n_eph"), bp::arg("multi_objective"),
                  bp::arg("parallel_batches") = 0u)),
             symbolic_regression_init_doc().c_str());

    // Making data from the gym available in python
    expose_data_from_the_gym<&gym::generate_koza_quintic>("generate_koza_quintic", generate_koza_quintic_doc());
    // From Our paper
    expose_data_from_the_gym<&gym::generate_P1>("generate_P1", generate_P1_doc());
    expose_data_from_the_gym<&gym::generate_P2>("generate_P2", generate_P2_doc());
    expose_data_from_the_gym<&gym::generate_P3>("generate_P3", generate_P3_doc());
    expose_data_from_the_gym<&gym::generate_P4>("generate_P4", generate_P4_doc());
    expose_data_from_the_gym<&gym::generate_P5>("generate_P5", generate_P5_doc());
    expose_data_from_the_gym<&gym::generate_P6>("generate_P6", generate_P6_doc());
    expose_data_from_the_gym<&gym::generate_P7>("generate_P7", generate_P7_doc());
    // From Vladi paper
    expose_data_from_the_gym<&gym::generate_kotanchek>("generate_kotanchek", generate_kotanchek_doc());
    expose_data_from_the_gym<&gym::generate_salutowicz>("generate_salutowicz", generate_salutowicz_doc());
    expose_data_from_the_gym<&gym::generate_salutowicz2d>("generate_salutowicz2d", generate_salutowicz2d_doc());
    expose_data_from_the_gym<&gym::generate_uball5d>("generate_uball5d", generate_uball5d_doc());
    expose_data_from_the_gym<&gym::generate_ratpol3d>("generate_ratpol3d", generate_ratpol3d_doc());
    expose_data_from_the_gym<&gym::generate_sinecosine>("generate_sinecosine", generate_sinecosine_doc());
    expose_data_from_the_gym<&gym::generate_ripple>("generate_ripple", generate_ripple_doc());
    expose_data_from_the_gym<&gym::generate_ratpol2d>("generate_ratpol2d", generate_ratpol2d_doc());
    // NIST data
    expose_data_from_the_gym<&gym::generate_chwirut1>("generate_chwirut1", generate_chwirut1_doc());
    expose_data_from_the_gym<&gym::generate_chwirut2>("generate_chwirut2", generate_chwirut2_doc());
    expose_data_from_the_gym<&gym::generate_daniel_wood>("generate_daniel_wood", generate_daniel_wood_doc());
    expose_data_from_the_gym<&gym::generate_gauss1>("generate_gauss1", generate_gauss1_doc());
    expose_data_from_the_gym<&gym::generate_kirby2>("generate_kirby2", generate_kirby2_doc());
    expose_data_from_the_gym<&gym::generate_lanczos2>("generate_lanczos2", generate_lanczos2_doc());
    expose_data_from_the_gym<&gym::generate_misra1b>("generate_misra1b", generate_misra1b_doc());
}
} // namespace dcgpy