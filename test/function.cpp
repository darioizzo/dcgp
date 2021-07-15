#define BOOST_TEST_MODULE dcgp_function_test
#include <boost/test/included/unit_test.hpp>

#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/algorithm/string/predicate.hpp>

#include <pagmo/threading.hpp>

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
// Workaround for a possible bug in MSVC
#if !defined(_MSC_VER) || defined(__clang__)
    BOOST_CHECK(!f1.is<void (*)(int)>());
    BOOST_CHECK(static_cast<const function<void()> &>(f1).extract<void (*)(int)>() == nullptr);
#endif
    BOOST_CHECK(static_cast<const function<void()> &>(f1).extract<void (*)()>() != nullptr);
    BOOST_CHECK_EXCEPTION(f1(), std::runtime_error, [](const std::runtime_error &re) {
        return boost::contains(
            re.what(),
            "This dcgp::function object cannot be invoked because it contains a null pointer to a C++ function");
    });
    BOOST_CHECK(f1.get_thread_safety() == pagmo::thread_safety::basic);

    // Copy construction.
    auto f1_copy(f1);
    BOOST_CHECK(f1_copy.is_valid());
    BOOST_CHECK(f1_copy.is<void (*)()>());
#if !defined(_MSC_VER) || defined(__clang__)
    BOOST_CHECK(!f1_copy.is<void (*)(int)>());
    BOOST_CHECK(static_cast<const function<void()> &>(f1_copy).extract<void (*)(int)>() == nullptr);
#endif
    BOOST_CHECK(static_cast<const function<void()> &>(f1_copy).extract<void (*)()>() != nullptr);
    BOOST_CHECK_EXCEPTION(f1_copy(), std::runtime_error, [](const std::runtime_error &re) {
        return boost::contains(
            re.what(),
            "This dcgp::function object cannot be invoked because it contains a null pointer to a C++ function");
    });
    BOOST_CHECK(f1_copy.get_thread_safety() == pagmo::thread_safety::basic);

    // Move construction.
    auto f1_move(std::move(f1_copy));
    BOOST_CHECK(f1_move.is_valid());
    BOOST_CHECK(f1_move.is<void (*)()>());
#if !defined(_MSC_VER) || defined(__clang__)
    BOOST_CHECK(!f1_move.is<void (*)(int)>());
    BOOST_CHECK(static_cast<const function<void()> &>(f1_move).extract<void (*)(int)>() == nullptr);
#endif
    BOOST_CHECK(static_cast<const function<void()> &>(f1_move).extract<void (*)()>() != nullptr);
    BOOST_CHECK_EXCEPTION(f1_move(), std::runtime_error, [](const std::runtime_error &re) {
        return boost::contains(
            re.what(),
            "This dcgp::function object cannot be invoked because it contains a null pointer to a C++ function");
    });
    BOOST_CHECK(f1_move.get_thread_safety() == pagmo::thread_safety::basic);
    BOOST_CHECK(!f1_copy.is_valid());

    // Copy assignment.
    f1_copy = f1_move;
    BOOST_CHECK(f1_copy.is_valid());
    BOOST_CHECK(f1_copy.is<void (*)()>());
#if !defined(_MSC_VER) || defined(__clang__)
    BOOST_CHECK(!f1_copy.is<void (*)(int)>());
    BOOST_CHECK(static_cast<const function<void()> &>(f1_copy).extract<void (*)(int)>() == nullptr);
#endif
    BOOST_CHECK(static_cast<const function<void()> &>(f1_copy).extract<void (*)()>() != nullptr);
    BOOST_CHECK_EXCEPTION(f1_copy(), std::runtime_error, [](const std::runtime_error &re) {
        return boost::contains(
            re.what(),
            "This dcgp::function object cannot be invoked because it contains a null pointer to a C++ function");
    });
    BOOST_CHECK(f1_copy.get_thread_safety() == pagmo::thread_safety::basic);

    // Move assignment.
    f1_move = std::move(f1_copy);
    BOOST_CHECK(f1_move.is_valid());
    BOOST_CHECK(f1_move.is<void (*)()>());
#if !defined(_MSC_VER) || defined(__clang__)
    BOOST_CHECK(!f1_move.is<void (*)(int)>());
    BOOST_CHECK(static_cast<const function<void()> &>(f1_move).extract<void (*)(int)>() == nullptr);
#endif

    BOOST_CHECK(static_cast<const function<void()> &>(f1_move).extract<void (*)()>() != nullptr);
    BOOST_CHECK_EXCEPTION(f1_move(), std::runtime_error, [](const std::runtime_error &re) {
        return boost::contains(
            re.what(),
            "This dcgp::function object cannot be invoked because it contains a null pointer to a C++ function");
    });
    BOOST_CHECK(f1_move.get_thread_safety() == pagmo::thread_safety::basic);
    BOOST_CHECK(!f1_copy.is_valid());

    // The simplest function.
    f1 = function<void()>{hello_world};
    BOOST_CHECK(f1.is_valid());
    f1();

    // A couple of other simple functions.
    function<double(double, double)> f2(double_add);
    BOOST_CHECK(f2(1, 2) == 3);
    function<double(const double &, const double &)> f3(double_add_ref);
    BOOST_CHECK(f3(1, 2) == 3);

    // Try with a functor too, setting a custom thread safety.
    struct local_f {
        bool operator()() const
        {
            return true;
        }
        pagmo::thread_safety get_thread_safety() const
        {
            return pagmo::thread_safety::none;
        }
    };
    function<bool()> f4(local_f{});
    BOOST_CHECK(f4());
    BOOST_CHECK(f4.get_thread_safety() == pagmo::thread_safety::none);
    BOOST_CHECK(f4.is<local_f>());
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
