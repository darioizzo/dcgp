#define BOOST_TEST_MODULE dcgp_symbolic_regression_test
#include <boost/test/unit_test.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/algorithms/sga.hpp>
#include <pagmo/algorithms/gaco.hpp>


#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(symbolic_regression_set)
{

    symbolic_regression udp({{1., 2.}, {0.3, -0.32}}, {{3./2.}, {0.02/0.32}});
    pagmo::problem prob{udp};
    pagmo::print(prob);
    pagmo::population pop{prob, 100};
    pagmo::sga uda(1000);
    pagmo::algorithm algo{uda};
    pagmo::print(algo);
    algo.set_verbosity(1u);
    pop = algo.evolve(pop);
    pagmo::print("Final Error: ", pop.champion_f(), "\n");
    pagmo::print("Expression: ", udp.pretty(pop.champion_x()));
}
