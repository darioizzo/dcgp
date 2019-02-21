#define BOOST_TEST_MODULE dcgp_expression_ann_test
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <stdexcept>

#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>
using namespace dcgp;

BOOST_AUTO_TEST_CASE(pdiv)
{
    // Initialize kernel
    kernel_set<double> pdiv_set({"pdiv"});
    kernel<double> pdiv = pdiv_set[0];

    BOOST_CHECK_EQUAL(pdiv({0.4, 0.5}), 0.8);
    BOOST_CHECK_EQUAL(pdiv({0., 0.5}), 0.);
    BOOST_CHECK_EQUAL(pdiv({0.4, 0.}), 1.);
    BOOST_CHECK(std::isfinite(pdiv({1., 1.2e-38})));
    BOOST_CHECK(std::isfinite(pdiv({1.2e-38, 1.})));
}
