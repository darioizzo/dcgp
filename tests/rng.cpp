#define BOOST_TEST_MODULE dcgp_differentiation_test
#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <pagmo/s11n.hpp>
#include <thread>
#include <vector>

#include <dcgp/rng.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(set_seed_and_next)
{
    // We check that the first N pseudo random numbers are identical if generated
    // right after the same seed is set and different otherwise.

    // We choose two seeds
    unsigned seed{0u}, seed2{1u};

    // Length of the pseudo-random sequence tested
    unsigned N = 10000u;

    // We generate three pseudo random sequences, two with the same seed
    random_device::set_seed(seed);
    std::vector<detail::random_engine_type::result_type> prs1;
    std::generate_n(std::back_inserter(prs1), N, random_device::next);

    random_device::set_seed(seed);
    std::vector<detail::random_engine_type::result_type> prs2;
    std::generate_n(std::back_inserter(prs2), N, random_device::next);

    random_device::set_seed(seed2);
    std::vector<detail::random_engine_type::result_type> prs3;
    std::generate_n(std::back_inserter(prs3), N, random_device::next);

    // We check that prs1 and prs2 are equal, since the seed was the same
    BOOST_CHECK(std::equal(prs1.begin(), prs1.end(), prs2.begin()));
    // We check that prs1 are prs3 are different since the seed was different
    BOOST_CHECK(!std::equal(prs1.begin(), prs1.end(), prs3.begin()));
}

// This test just runs calls to random_device::next() in two separate threads. If this executable
// is compiled with -fsanitize=thread in clang/gcc, it should check that the locking logic
// in random_device is correct.
BOOST_AUTO_TEST_CASE(data_races_test)
{
    unsigned N = 10000u;
    std::vector<detail::random_engine_type::result_type> prs4, prs5;
    std::thread t1([&]() { std::generate_n(std::back_inserter(prs4), N, random_device::next); });
    std::thread t2([&]() { std::generate_n(std::back_inserter(prs5), N, random_device::next); });
    t1.join();
    t2.join();
}

BOOST_AUTO_TEST_CASE(rng_serialization_test)
{
    const int ntrials = 100;
    std::mt19937 rng;
    using r_type = std::mt19937;
    using ia_type = boost::archive::binary_iarchive;
    using oa_type = boost::archive::binary_oarchive;
    auto rng_save = [](const r_type &r) {
        std::stringstream ss;
        {
            oa_type oarchive(ss);
            oarchive << r;
        }
        return ss.str();
    };
    auto rng_load = [](const std::string &str, r_type &r) {
        std::stringstream ss;
        ss.str(str);
        {
            ia_type iarchive(ss);
            iarchive >> r;
        }
    };
    std::uniform_int_distribution<r_type::result_type> dist;
    for (auto i = 0; i < ntrials; ++i) {
        auto seed = dist(rng);
        r_type r;
        r.seed(seed);
        auto str = rng_save(r);
        std::vector<r_type::result_type> v1;
        std::generate_n(std::back_inserter(v1), 100, r);
        auto r_copy(r);
        rng_load(str, r);
        std::vector<r_type::result_type> v2;
        std::generate_n(std::back_inserter(v2), 100, r);
        BOOST_CHECK_EQUAL_COLLECTIONS(v1.begin(), v1.end(), v2.begin(), v2.end());
        BOOST_CHECK(r_copy == r);
    }
}
