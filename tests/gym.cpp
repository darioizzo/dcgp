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
    gym::generate_vladi1(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi1(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi2(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi2(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi3(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi3(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi4(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi4(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi5(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi5(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi6(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi6(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi7(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi7(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
    gym::generate_vladi8(points, labels);
    for (decltype(points.size()) i = 0u; i < points.size(); ++i) {

        BOOST_CHECK_EQUAL(labels[i][0], gym::detail::vladi8(points[i]));
        BOOST_CHECK_EQUAL(labels.size(), points.size());
    }
}
