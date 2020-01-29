#define BOOST_TEST_MODULE dcgp_expression_weighted_test
#include <boost/test/included/unit_test.hpp>

#include <dcgp/expression_weighted.hpp>
#include <dcgp/kernel_set.hpp>

#include "helpers.hpp"

using namespace dcgp;


BOOST_AUTO_TEST_CASE(get_set_weight_test)
{
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::vector<unsigned> arity{3, 2, 1};
    expression_weighted<double> ex(1, 1, 2, 3, 4, arity, basic_set(), 0u);
    // we use the float of the input_id as weight
    std::vector<double> weights{0., 1., 2., 0., 1., 2., 0., 1., 0., 1., 0., 0.};
    ex.set_weights(weights);
    // check getting single weights
    BOOST_CHECK(ex.get_weight(1u, 0u) == 0.);
    BOOST_CHECK(ex.get_weight(1u, 1u) == 1.);
    BOOST_CHECK(ex.get_weight(1u, 2u) == 2.);
    // skip check for node 2 and 3, 2 is same as 1; 3 is same as 4
    BOOST_CHECK(ex.get_weight(4u, 0u) == 0.);
    BOOST_CHECK(ex.get_weight(4u, 1u) == 1.);
    // check node 5 and 6 (they are the same)
    BOOST_CHECK(ex.get_weight(5u, 0u) == 0.);
    BOOST_CHECK(ex.get_weight(6u, 0u) == 0.);
    // change single weights in expression and reference weights
    ex.set_weight(1u, 0u, 10.);
    weights[0] = 10.;
    ex.set_weight(4u, 1u, 10.);
    // position in reference weights is r*arity(c_1) + r*arity(c_2) = 10
    weights[9] = 10.;
    // check all weights together
    CHECK_CLOSE_V(ex.get_weights(), weights, 1e-12);
}