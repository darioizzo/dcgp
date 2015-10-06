#include <vector>
#define BOOST_TEST_MODULE dcgp_differentiation_test
#include <boost/test/unit_test.hpp>

#include "../src/expression.hpp"
#include "../src/function_set.hpp"

using namespace dcgp;

BOOST_AUTO_TEST_CASE(differentiation_basic_set)
{
    function_set basic_set({"sum","diff","mul","div"});
    expression ex(3, 1, 1, 20, 21, 2, basic_set(), 0);

    std::vector<double> in({1, 1, 1});
    std::vector<double> in2({-1, 1, 1});

    // We set the expression to 2 y^2 (x + z)^2
    ex.set({0, 2, 1, 0, 1, 3, 0, 2, 0, 2, 5, 1, 0, 6, 6, 1, 3, 7, 2, 5, 3, 3, 0, 3, 0, 10, 10, 2, 1, 7, 3, 0, 2, 2, 12, 5, 2, 0, 11, 0, 5, 5, 1, 4, 7, 2, 0, 13, 1, 2, 12, 0, 2, 6, 3, 20, 6, 1, 13, 11, 14});
    // We compute the jet of derivatives up to order 2
    auto jet = ex.taylor(in, 2);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({0, 0, 0}), 8);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({1, 0, 0}), 8);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({0, 1, 0}), 16);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({0, 0, 1}), 8);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({2, 0, 0}), 4);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({0, 2, 0}), 16);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({0, 0, 2}), 4);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({1, 1, 0}), 16);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({0, 1, 1}), 16);
    BOOST_CHECK_EQUAL(jet[0].get_derivative({1, 0, 1}), 4);

    // We set the expression to y / (x^2 z(x + y))
    ex.set({3, 1, 0, 0, 1, 0, 0, 4, 4, 3, 3, 4, 2, 0, 2, 2, 7, 6, 1, 7, 2, 0, 0, 5, 0, 9, 9, 3, 6, 4, 2, 12, 1, 2, 10, 10, 0, 10, 12, 0, 14, 13, 1, 3, 14, 3, 6, 7, 2, 15, 1, 0, 19, 13, 1, 8, 3, 3, 2, 18, 18});
    auto jet1 = ex.taylor(in, 2);

    BOOST_CHECK_EQUAL(jet1[0].get_derivative({0, 0, 0}), 0.5);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({1, 0, 0}), - 1.25);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({0, 1, 0}), 0.25);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({0, 0, 1}), -0.5);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({2, 0, 0}), 4.25);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({0, 2, 0}), -0.25);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({0, 0, 2}), 1);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({1, 1, 0}), -0.5);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({0, 1, 1}), -0.25);
    BOOST_CHECK_EQUAL(jet1[0].get_derivative({1, 0, 1}), 1.25);

    // We set the expression to 2 / (y (z-x) )
    ex.set({0, 2, 2, 1, 2, 0, 3, 3, 4, 3, 5, 1, 3, 6, 2, 0, 4, 1, 1, 8, 2, 2, 8, 2, 3, 0, 7, 3, 4, 2, 0, 2, 1, 0, 3, 8, 0, 0, 3, 3, 12, 11, 3, 10, 4, 3, 4, 5, 0, 6, 5, 0, 8, 12, 1, 11, 6, 0, 4, 15, 7});
    auto jet2 = ex.taylor(in2, 6);
    BOOST_CHECK_EQUAL(jet2[0].get_derivative({0, 0, 0}), 1);
    BOOST_CHECK_EQUAL(jet2[0].get_derivative({1, 2, 3}), - 3);
    BOOST_CHECK_EQUAL(jet2[0].get_derivative({3, 2, 1}), - 3);
    BOOST_CHECK_EQUAL(jet2[0].get_derivative({2, 3, 1}), 4.5);
}