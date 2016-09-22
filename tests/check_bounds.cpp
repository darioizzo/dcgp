#include <random>
#define BOOST_TEST_MODULE dcgp_check_bounds_test
#include <boost/test/unit_test.hpp>

#include "../include/dcgp.hpp"
#include "helpers.hpp"

using namespace dcgp;

BOOST_AUTO_TEST_CASE(check_bounds)
{
    // Random seed
    std::random_device rd;

    function_set basic_set({"sum","diff","mul","div"});

    // Testing an expression with arity 3, levels-back 2
    expression ex(3,1,2,3,2,3,basic_set(),rd());

    CHECK_EQUAL_V(ex.get_lb(), std::vector<unsigned int>({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,0,3,3,3,5}));
    CHECK_EQUAL_V(ex.get_ub(), std::vector<unsigned int>({3,2,2,2,3,2,2,2,3,4,4,4,3,4,4,4,3,6,6,6,3,6,6,6,8}));
}
