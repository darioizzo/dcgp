#define BOOST_TEST_MODULE dcgp_wrapped_functions_test
#include <boost/test/included/unit_test.hpp>

#include <cmath>
#include <initializer_list>
#include <sstream>
#include <vector>

#include <dcgp/function.hpp>
#include <dcgp/s11n.hpp>
#include <dcgp/wrapped_functions.hpp>
#include <dcgp/wrapped_functions_s11n_implement.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(pdiv)
{

    // test with arity of 2
    {
        std::vector<double> v({0.4, 0.5});

        BOOST_CHECK_EQUAL(my_pdiv<double>(v), 0.8);

        v[0] = 0;
        v[1] = 0.5;
        BOOST_CHECK_EQUAL(my_pdiv<double>(v), 0.);

        v[0] = 0.4;
        v[1] = 0.;
        BOOST_CHECK_EQUAL(my_pdiv<double>(v), 1.);

        v[0] = 1.;
        v[1] = 1.2e-38;
        BOOST_CHECK(std::isfinite(my_pdiv<double>(v)));

        v[0] = 1.2e-38;
        v[1] = 1.;
        BOOST_CHECK(std::isfinite(my_pdiv<double>(v)));
    }
    // test with arity 5
    {
        std::vector<double> v({0.5, 0.4, 0.2, 0.2, 0.1});

        BOOST_CHECK_CLOSE(my_pdiv<double>(v), 312.5, 1e-12);

        v[3] = 0.;
        BOOST_CHECK_EQUAL(my_pdiv<double>(v), 1.);

        v[3] = 1.2e-38;
        BOOST_CHECK(std::isfinite(my_pdiv<double>(v)));
    }
    // Serialization test.
    {
        function<double(const std::vector<double> &)> f{my_pdiv<double>};
        std::stringstream ss;
        {
            boost::archive::binary_oarchive oarchive(ss);
            oarchive << f;
        }
        f = function<double(const std::vector<double> &)>{};
        BOOST_CHECK(!f.is<my_pdiv_func<double>>());
        {
            boost::archive::binary_iarchive iarchive(ss);
            iarchive >> f;
        }
        BOOST_CHECK(f.is<my_pdiv_func<double>>());
        BOOST_CHECK(f({2., 3.}) == my_pdiv<double>({2., 3.}));
    }
}

BOOST_AUTO_TEST_CASE(my_sqrt_test)
{
    // test with arity of 2
    {
        std::vector<double> v({4, 0.34});

        BOOST_CHECK_CLOSE(my_sqrt<double>(v), 2, 1e-4);

        v[0] = 0;
        v[1] = 0.5;
        BOOST_CHECK_CLOSE(my_sqrt<double>(v), 0., 1e-4);

        v[0] = 16;
        v[1] = 0.;
        BOOST_CHECK_CLOSE(my_sqrt<double>(v), 4., 1e-4);

        v[0] = 1.;
        v[1] = 1.2e-38;
        BOOST_CHECK(std::isfinite(my_sqrt<double>(v)));

        v[0] = 1.2e-38;
        v[1] = 1.;
        BOOST_CHECK(std::isfinite(my_sqrt<double>(v)));

        // This fails in MinGW 6.2 as a known bug So for the time being we deactivate it as
        // our appveyor build are using that MinGW
        // v[0] = -1.2e-38;
        // v[1] = 1.;
        // std::cout << "Res: " << my_sqrt(v) << std::endl;
        // BOOST_CHECK(!std::isfinite(my_sqrt(v)));
    }
    // test with arity 5
    {
        std::vector<double> v{4, 0.4, 0.2, 0.2, 0.1};

        BOOST_CHECK_CLOSE(my_sqrt<double>(v), 2., 1e-4);

        v[3] = 0.;
        BOOST_CHECK_CLOSE(my_sqrt<double>(v), 2., 1e-4);

        v[3] = 1.2e-38;
        BOOST_CHECK(std::isfinite(my_sqrt<double>(v)));
    }
    // Serialization test.
    {
        function<double(const std::vector<double> &)> f{my_sqrt<double>};
        std::stringstream ss;
        {
            boost::archive::binary_oarchive oarchive(ss);
            oarchive << f;
        }
        f = function<double(const std::vector<double> &)>{};
        BOOST_CHECK(!f.is<my_sqrt_func<double>>());
        {
            boost::archive::binary_iarchive iarchive(ss);
            iarchive >> f;
        }
        BOOST_CHECK(f.is<my_sqrt_func<double>>());
        BOOST_CHECK(f({2.}) == std::sqrt(2.));
    }
}

BOOST_AUTO_TEST_CASE(my_gaussian_test)
{

    // test with arity of 2
    {
        std::vector<double> v({0.132, 0.34});
        BOOST_CHECK_CLOSE(my_gaussian<double>(v), 0.98272692007, 1e-4);
    }
    // test with arity 5
    {
        std::vector<double> v({0.132, 0.4, 0.2, 0.2, 0.1});
        BOOST_CHECK_CLOSE(my_gaussian<double>(v), 0.98272692007, 1e-4);
    }
    // Serialization test.
    {
        function<double(const std::vector<double> &)> f{my_gaussian<double>};
        std::stringstream ss;
        {
            boost::archive::binary_oarchive oarchive(ss);
            oarchive << f;
        }
        f = function<double(const std::vector<double> &)>{};
        BOOST_CHECK(!f.is<my_gaussian_func<double>>());
        {
            boost::archive::binary_iarchive iarchive(ss);
            iarchive >> f;
        }
        BOOST_CHECK(f.is<my_gaussian_func<double>>());
        BOOST_CHECK(f({0.132, 0.4, 0.2, 0.2, 0.1}) == my_gaussian<double>({0.132, 0.4, 0.2, 0.2, 0.1}));
    }
}
