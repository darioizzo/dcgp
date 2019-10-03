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

