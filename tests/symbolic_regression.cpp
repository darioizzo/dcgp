#define BOOST_TEST_MODULE dcgp_symbolic_regression_test
#include <boost/test/unit_test.hpp>

#include <dcgp/problems/symbolic_regression.hpp>
#include <dcgp/kernel_set.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(symbolic_regression_set)
{
    symbolic_regression udp{};
    BOOST_CHECK(udp.fitness({1.,2.}) == pagmo::vector_double{1.});

}
