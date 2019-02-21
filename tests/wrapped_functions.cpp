#define BOOST_TEST_MODULE dcgp_wrapped_functions_test
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <dcgp/wrapped_functions.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(pdiv)
{

    // test with arity of 2
    {
        std::vector<double> v({0.4, 0.5});

        BOOST_CHECK_EQUAL(my_pdiv(v), 0.8);

        v[0] = 0;
        v[1] = 0.5;
        BOOST_CHECK_EQUAL(my_pdiv(v), 0.);

        v[0] = 0.4;
        v[1] = 0.;
        BOOST_CHECK_EQUAL(my_pdiv(v), 1.);

        v[0] = 1.;
        v[1] = 1.2e-38;
        BOOST_CHECK(std::isfinite(my_pdiv(v)));

        v[0] = 1.2e-38;
        v[1] = 1.;
        BOOST_CHECK(std::isfinite(my_pdiv(v)));
    }
    // test with arity 5
    {
        std::vector<double> v({0.5, 0.4, 0.2, 0.2, 0.1});

        BOOST_CHECK_CLOSE(my_pdiv(v), 312.5, 1e-12);

        v[3] = 0.;
        BOOST_CHECK_EQUAL(my_pdiv(v), 1.);

        v[3] = 1.2e-38;
        BOOST_CHECK(std::isfinite(my_pdiv(v)));
    }
}
