#define BOOST_TEST_MODULE dcgp_function_test
#include <boost/test/unit_test.hpp>

#include <sstream>

#include <dcgp/function.hpp>
#include <dcgp/s11n.hpp>

using namespace dcgp;

struct hello_world_func {
    void operator()() const {}
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

inline constexpr auto hello_world = hello_world_func{};

DCGP_S11N_FUNCTION_EXPORT(hwf, hello_world_func, void)

BOOST_AUTO_TEST_CASE(basic_test)
{
    function<void()> f{hello_world};
    f();
}

BOOST_AUTO_TEST_CASE(function_serialization_test)
{
    function<void()> f{hello_world};
    std::stringstream ss;
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << f;
    }
    f = function<void()>{};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> f;
    }
    f();
}
