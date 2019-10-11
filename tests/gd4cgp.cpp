#define BOOST_TEST_MODULE dcgp_es4cgp_test
#include <boost/test/unit_test.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/rosenbrock.hpp>

#include <dcgp/algorithms/gd4cgp.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction_test)
{
    BOOST_CHECK_NO_THROW(gd4cgp(1u, 0.1, 0.001));
    BOOST_CHECK_THROW(gd4cgp(1u, 0.1, 0.2), std::invalid_argument);
    BOOST_CHECK_THROW(gd4cgp(1u, -0.1, 0.001), std::invalid_argument);
    BOOST_CHECK_THROW(gd4cgp(1u, 0.1, -0.001), std::invalid_argument);
    BOOST_CHECK_THROW(gd4cgp(1u, -0.1, -1.001), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(evolve_test)
{
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    symbolic_regression udp({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, 1, 10, 3, 2, basic_set(), 2u, 0u);
    pagmo::problem prob{udp};
    pagmo::population pop1{prob, 5u};
    gd4cgp uda1{10, 1., 1e-3};
    uda1.set_verbosity(1u);
    pop1 = uda1.evolve(pop1);
}

