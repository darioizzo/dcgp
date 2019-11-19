#define BOOST_TEST_MODULE dcgp_es4cgp_test
#include <boost/test/included/unit_test.hpp>

#include <vector>

#include <dcgp/gym.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(gym_test)
{
    std::vector<std::vector<double>> points, labels;
    gym::generate_koza_quintic(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {
        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::koza_quintic(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P1(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P1(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P2(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P2(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P3(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P3(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P4(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P4(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P5(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P5(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P6(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P6(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_P7(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::P7(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_kotanchek(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::kotanchek(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_salutowicz(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::salutowicz(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_salutowicz2d(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::salutowicz2d(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_uball5d(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::uball5d(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_ratpol3d(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::ratpol3d(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_sinecosine(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::sinecosine(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_ripple(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::ripple(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_ratpol2d(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::ratpol2d(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
}
