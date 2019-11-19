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

BOOST_AUTO_TEST_CASE(construction_test)
{
// We test that kernels can be constructed with the shipped functions
BOOST_CHECK_NO_THROW(kernel<double>(my_sum<double>, print_my_sum, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_diff<double>, print_my_diff, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_mul<double>, print_my_mul, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_div<double>, print_my_div, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_pdiv<double>, print_my_pdiv, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_sin<double>, print_my_sin, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_cos<double>, print_my_cos, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_tanh<double>, print_my_tanh, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_exp<double>, print_my_exp, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_gaussian<double>, print_my_gaussian, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_sqrt<double>, print_my_sqrt, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_psqrt<double>, print_my_sqrt, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_log<double>, print_my_log, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_sig<double>, print_my_sig, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_relu<double>, print_my_relu, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_elu<double>, print_my_elu, "w"));
BOOST_CHECK_NO_THROW(kernel<double>(my_isru<double>, print_my_isru, "w"));
}

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
