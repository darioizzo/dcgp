#define BOOST_TEST_MODULE dcgp_es4cgp_test
#include <boost/test/unit_test.hpp>

#include <vector>

#include <dcgp/gym.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(gym_test)
{
    std::vector<std::vector<double>> points, labels;
    gym::generate_koza_quintic(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {
        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::koza_quintic_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P1(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P1_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P2(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P2_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P3(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P3_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P4(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P4_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P5(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P5_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P6(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P6_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P7(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P7_fun(points[i][0]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
}
