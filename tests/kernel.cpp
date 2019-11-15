#define BOOST_TEST_MODULE dcgp_kernel_test
#include <boost/test/included/unit_test.hpp>

#include <initializer_list>
#include <sstream>

#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/s11n.hpp>
#include <dcgp/wrapped_functions.hpp>
#include <dcgp/wrapped_functions_s11n_implement.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(s11n_test)
{
    {
        kernel<double> k1(my_sqrt<double>, print_my_sqrt, "carogna");
        kernel<double> k2(my_sin<double>, print_my_sin, "putrida");
        std::stringstream ss;
        {
            boost::archive::binary_oarchive oarchive(ss);
            oarchive << k1;
        }
        k1 = k2;
        BOOST_CHECK(k1.get_name() == "putrida");
        {
            boost::archive::binary_iarchive iarchive(ss);
            iarchive >> k1;
        }
        BOOST_CHECK(k1.get_name() == "carogna");
    }
    {
        kernel_set<double> ks1{{"diff", "mul", "div"}};
        kernel_set<double> ks2{{"diff"}};
        std::stringstream ss;
        {
            boost::archive::binary_oarchive oarchive(ss);
            oarchive << ks1;
        }
        ks1 = ks2;
        BOOST_CHECK(ks1().size() == 1u);
        {
            boost::archive::binary_iarchive iarchive(ss);
            iarchive >> ks1;
        }
        BOOST_CHECK(ks1().size() == 3u);
    }
}
