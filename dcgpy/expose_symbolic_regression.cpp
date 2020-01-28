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

#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/algorithms/gd4cgp.hpp>
#include <dcgp/algorithms/mes4cgp.hpp>
#include <dcgp/algorithms/momes4cgp.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/kernel.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

#include <pagmo/algorithm.hpp>
#include <pagmo/problem.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"

PAGMO_S11N_ALGORITHM_IMPLEMENT(dcgp::es4cgp)
PAGMO_S11N_ALGORITHM_IMPLEMENT(dcgp::mes4cgp)
PAGMO_S11N_ALGORITHM_IMPLEMENT(dcgp::momes4cgp)
PAGMO_S11N_ALGORITHM_IMPLEMENT(dcgp::gd4cgp)
PAGMO_S11N_PROBLEM_IMPLEMENT(dcgp::symbolic_regression)

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
                     bool multi_objective, unsigned parallel_batches, const std::string loss_s) {
                     auto kernels_v = l_to_v<kernel<double>>(kernels);
                     auto vvd_points = to_vv<double>(points);
                     auto vvd_labels = to_vv<double>(labels);
                     return ::new dcgp::symbolic_regression(vvd_points, vvd_labels, rows, cols, levels_back, arity,
                                                            kernels_v, n_eph, multi_objective, parallel_batches, loss_s);
                 },
                 bp::default_call_policies(),
                 (bp::arg("points"), bp::arg("labels"), bp::arg("rows") = 1, bp::arg("cols") = 16,
                  bp::arg("levels_back") = 17, bp::arg("arity") = 2, bp::arg("kernels"), bp::arg("n_eph") = 0u,
                  bp::arg("multi_objective") = false, bp::arg("parallel_batches") = 0u, bp::arg("loss") = "MSE")),
             symbolic_regression_init_doc().c_str())
        .def(
            "pretty", +[](const dcgp::symbolic_regression &instance,
                          const bp::object &x) { return instance.pretty(to_v<double>(x)); })
        .def(
            "prettier", +[](const dcgp::symbolic_regression &instance,
                            const bp::object &x) { return instance.prettier(to_v<double>(x)); })
        .def(
            "predict",
            +[](const dcgp::symbolic_regression &instance, const bp::object &points, const bp::object &x) {
                try { //
                    return vvector_to_ndarr(instance.predict(to_vv<double>(points), to_v<double>(x)));
                } catch (...) {
                    PyErr_Clear();
                    return vector_to_ndarr(instance.predict(to_v<double>(points), to_v<double>(x)));
                }
            },
            (bp::arg("points"), bp::arg("chromosome")), symbolic_regression_predict_doc().c_str())
        .def("__repr__", &dcgp::symbolic_regression::get_extra_info);

    // We expose the UDAs
    // ES-4CGP (Evolutionary Strategy for Cartesian Genetic Programming)
    auto es4cgp_ = pg::expose_algorithm<dcgp::es4cgp>("es4cgp", es4cgp_doc().c_str());
    es4cgp_.def(bp::init<unsigned, unsigned, double, bool>(
        (bp::arg("gen") = 1u, bp::arg("mut_n") = 1u, bp::arg("ftol") = 1e-4, bp::arg("learn_constants") = true)));
    es4cgp_.def(bp::init<unsigned, unsigned, double, bool, unsigned>(
        (bp::arg("gen") = 1u, bp::arg("mut_n") = 1u, bp::arg("ftol") = 1e-4, bp::arg("learn_constants") = true,
         bp::arg("seed"))));
    es4cgp_.def("get_seed", &es4cgp::get_seed, generic_uda_get_seed_doc().c_str());
    // es4cgp_ needs an ad hoc exposition for the log as one entry is a vector (constants)
    es4cgp_.def(
        "get_log",
        +[](const dcgp::es4cgp &a) -> bp::list {
            bp::list retval;
            for (const auto &t : a.get_log()) {
                retval.append(bp::make_tuple(std::get<0>(t), std::get<1>(t), std::get<2>(t),
                                             pygmo::vector_to_ndarr(std::get<3>(t)), std::get<4>(t)));
            }
            return retval;
        },
        es4cgp_get_log_doc().c_str());
    // GD-4CGP (Gradient Descent for Cartesian Genetic Programming)
    auto gd4cgp_ = pg::expose_algorithm<dcgp::gd4cgp>("gd4cgp", gd4cgp_doc().c_str());
    gd4cgp_.def(
        bp::init<unsigned, double, double>((bp::arg("max_iter") = 1u, bp::arg("lr") = 1., bp::arg("lr_min") = 1e-3)));
    pg::expose_algo_log(gd4cgp_, gd4cgp_get_log_doc().c_str());
    // MES-4CGP (Memetic Evolutionary Strategy for Cartesian Genetic Programming)
    auto mes4cgp_ = pg::expose_algorithm<dcgp::mes4cgp>("mes4cgp", mes4cgp_doc().c_str());
    mes4cgp_.def(
        bp::init<unsigned, unsigned, double>((bp::arg("gen") = 1u, bp::arg("mut_n") = 1u, bp::arg("ftol") = 1e-4)));
    mes4cgp_.def(bp::init<unsigned, unsigned, double, unsigned>(
        (bp::arg("gen") = 1u, bp::arg("mut_n") = 1u, bp::arg("ftol") = 1e-4, bp::arg("seed"))));
    mes4cgp_.def("get_seed", &dcgp::mes4cgp::get_seed, generic_uda_get_seed_doc().c_str());
    // mes4cgp_ needs an ad hoc exposition for the log as one entry is a vector (constants)
    mes4cgp_.def(
        "get_log",
        +[](const dcgp::mes4cgp &a) -> bp::list {
            bp::list retval;
            for (const auto &t : a.get_log()) {
                retval.append(bp::make_tuple(std::get<0>(t), std::get<1>(t), std::get<2>(t),
                                             pygmo::vector_to_ndarr(std::get<3>(t)), std::get<4>(t)));
            }
            return retval;
        },
        mes4cgp_get_log_doc().c_str());

    // MOMES-4CGP (Multi-Objective Memetic Evolutionary Strategy for Cartesian Genetic Programming)
    auto momes4cgp_ = pg::expose_algorithm<dcgp::momes4cgp>("momes4cgp", momes4cgp_doc().c_str());
    momes4cgp_.def(bp::init<unsigned, unsigned>((bp::arg("gen") = 1u, bp::arg("max_mut") = 1u)));
    momes4cgp_.def(
        bp::init<unsigned, unsigned, unsigned>((bp::arg("gen") = 1u, bp::arg("max_mut") = 1u, bp::arg("seed"))));
    momes4cgp_.def("get_seed", &dcgp::momes4cgp::get_seed, generic_uda_get_seed_doc().c_str());
    pg::expose_algo_log(momes4cgp_, momes4cgp_get_log_doc().c_str());

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