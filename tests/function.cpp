#define BOOST_TEST_MODULE dcgp_function_test
#include <boost/test/unit_test.hpp>

#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string/predicate.hpp>

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

double double_add(double a, double b)
{
    return a + b;
}

double double_add_ref(const double &a, const double &b)
{
    return a + b;
}

BOOST_AUTO_TEST_CASE(function_basic_tests)
{
    // Default construction.
    function<void()> f1;
    BOOST_CHECK(f1.is_valid());
    BOOST_CHECK(f1.is<void (*)()>());
    BOOST_CHECK(!f1.is<void (*)(int)>());
    BOOST_CHECK(static_cast<const function<void()> &>(f1).extract<void (*)()>() != nullptr);
    BOOST_CHECK(static_cast<const function<void()> &>(f1).extract<void (*)(int)>() == nullptr);
    BOOST_CHECK_EXCEPTION(f1(), std::runtime_error, [](const std::runtime_error &re) {
        return boost::contains(
            re.what(),
            "This dcp::function object cannot be invoked because it contains a null pointer to a C++ function");
    });

    // The simplest function.
    f1 = function<void()>{hello_world};
    BOOST_CHECK(f1.is_valid());
    f1();

    function<double(double, double)> f2(double_add);
    BOOST_CHECK(f2(1, 2) == 3);

    function<double(const double &, const double &)> f3(double_add_ref);
    BOOST_CHECK(f3(1, 2) == 3);
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
