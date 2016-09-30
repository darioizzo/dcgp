#include <vector>
#include <random>
#include <algorithm>
#define BOOST_TEST_MODULE dcgp_compute_test
#include <boost/test/unit_test.hpp>

#include "../include/expression.hpp"
#include "../include/function_set.hpp"

using namespace dcgp;

BOOST_AUTO_TEST_CASE(mutate)
{
    // Random seed
    std::random_device rd;
    // We define a d-CGP expression and get the reference chromosome
    function_set<double> basic_set({"sum","diff","mul","div"});
    // Number of random trials
    unsigned int N = 100;
    // A d-CGP expression
    expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());

    // We test mutate(idx). Does the idx gene change? Do all other stay the same?
    for (auto i = 0u; i < N; ++i)
    {
        std::vector<unsigned int> x = ex.get();
        auto idx = rd() % x.size();
        ex.mutate(idx);
        for (auto i = 0u; i < x.size(); ++i) {
            if (i == idx) {
                BOOST_CHECK(x[i] != ex.get()[i]);
            }
            else {
                BOOST_CHECK_EQUAL(x[i], ex.get()[i]);
            }
        }
    }

    // We test mutate_active. Was the mutated gene active?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i)
        {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_active();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned int idx = it.first - x.begin();
            auto ag = ex.get_active_genes();
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
        }
    }

    // We test mutate_active_fgene. Was the mutated gene active? Was it a function gene?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i)
        {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_active_fgene();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned int idx = it.first - x.begin();
            auto ag = ex.get_active_genes();
            // was it an active gene?
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
            // was it a function gene?
            BOOST_CHECK(idx < x.size() - ex.get_m());
            BOOST_CHECK_EQUAL(idx % 3, 0);
        }
    }

    // We test mutate_active_cgene. Was the mutated gene active? Was it a connection gene?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i)
        {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_active_cgene();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned int idx = it.first - x.begin();
            auto ag = ex.get_active_genes();
            // was it an active gene?
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
            // was it a connection gene?
            BOOST_CHECK(idx < x.size() - ex.get_m());
            BOOST_CHECK(((idx % 3) == 1 || (idx  % 3) == 2) == true);
        }
    }

    // We test mutate_ogene. Was the mutated gene active? Was it an output gene?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i)
        {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_ogene();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned int idx = it.first - x.begin();
            auto ag = ex.get_active_genes();
            // was it an active gene?
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
            // was it an output gene?
            BOOST_CHECK(idx >= x.size() - ex.get_m());
        }
    }
}
