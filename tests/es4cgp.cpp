#define BOOST_TEST_MODULE dcgp_es4cgp_test
#include <boost/test/unit_test.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(dcgp_es4cgp_test)
{
    symbolic_regression udp({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}});
    pagmo::problem prob{udp};
    pagmo::print(prob);
    pagmo::population pop{prob, 10};
    dcgp::es4cgp uda(100);
    pagmo::algorithm algo{uda};
    pagmo::print(algo);
    algo.set_verbosity(1u);
    pop = algo.evolve(pop);
}
