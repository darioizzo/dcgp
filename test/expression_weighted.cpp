#define BOOST_TEST_MODULE dcgp_expression_weighted_test
#include <boost/test/included/unit_test.hpp>

#include <dcgp/expression_weighted.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/wrapped_functions_s11n_implement.hpp>
#include <pagmo/s11n.hpp>
#include <pagmo/io.hpp>


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

BOOST_AUTO_TEST_CASE(s11n_test)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    expression_weighted<double> ex(2u, 2u, 20u, 20u, 3u, 2u, basic_set(), rd());
    // We change the weight values to test that non default values are deserialized
    ex.set_weights(std::vector<double>(ex.get_weights().size(), 0.123));
    const auto orig = boost::lexical_cast<std::string>(ex);

    std::stringstream ss;
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << ex;
    }
    ex = expression_weighted<double>(2, 2, 2, 2, 3, 2, basic_set(), rd());
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> ex;
    }

    BOOST_CHECK(orig == boost::lexical_cast<std::string>(ex));
}