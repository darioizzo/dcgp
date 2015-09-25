#include <iostream>
#include <iomanip>
#define BOOST_TEST_MODULE dcgp_compute_test
#include <boost/test/unit_test.hpp>

#include "../src/dcgp.h"
#include "helpers.hpp"

using namespace dcgp;

BOOST_AUTO_TEST_CASE(compute)
{
    function_set basic_set({"sum","diff","mul","div"});

    /// Testing over Miller's test case from the PPSN 2014 tutorial
    expression ex1(2,4,2,3,4,basic_set());
    ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});
    BOOST_CHECK_EQUAL(ex1({2.,3.}),  std::vector<double>({5,6,-18,0}));
    BOOST_CHECK_EQUAL(ex1({1.,-1.}), std::vector<double>({0,-1,-1,0}));
    BOOST_CHECK(EPSILON_COMPARE(ex1({-.123,2.345}), std::vector<double>({2.222,-0.288435,0.676380075,0}), 1e-14) == true);

    ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 6, 5, 7, 3});
    BOOST_CHECK_EQUAL(ex1({2.,-3.}), std::vector<double>({6,-6,-18,0}));
    BOOST_CHECK_EQUAL(ex1({1.,-1.}), std::vector<double>({2,-1,-1,0}));
    BOOST_CHECK(EPSILON_COMPARE(ex1({-.123,2.345}), std::vector<double>({-4.69,-0.288435,0.676380075,0}), 1e-14) == true);

    /// Testing over a single row program
    dcgp::expression ex2(4,1,1,10,10,basic_set());
    ex2.set({2, 3, 0, 0, 2, 2, 3, 0, 1, 1, 5, 4, 2, 6, 1, 0, 7, 7, 3, 6, 7, 1, 7, 6, 2, 4, 10, 2, 3, 2, 10});  ///(x/y)/(2z-(t*x))
    BOOST_CHECK(EPSILON_COMPARE(ex2({2.,3.,4.,-2.}), std::vector<double>({0.055555555555555552}), 1e-14) == true);
    BOOST_CHECK_EQUAL(ex2({-1.,1.,-1.,1.}), std::vector<double>({1}));
}
