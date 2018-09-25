#include <audi/audi.hpp>
#include <random>
#include <string>
#include <vector>

#define BOOST_TEST_MODULE dcgp_compute_test
#include <boost/test/unit_test.hpp>

#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>

#include "helpers.hpp"

using namespace dcgp;


double test_loss(unsigned int n, unsigned int m, unsigned int r, unsigned int c, unsigned int l, unsigned int a,
                 unsigned int N) // number of samples
{
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    dcgp::expression<double> ex(n, m, r, c, l, a, basic_set(), 123);

    // creates N data points
    std::default_random_engine re;
    std::vector<std::vector<double>> in;
    std::vector<std::vector<double>> out;
    std::vector<double> in_point(n);
    std::vector<double> out_point(m);
    for (auto i = 0u; i < N; ++i) {
        for (auto j = 0u; j < n; ++j) {
            in_point[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        }
        out_point = ex(in_point);
        in.push_back(in_point);
        out.push_back(out_point);
    }
    // computes the loss
    return ex.loss(in, out, "MSE", true);
}

audi::gdual_d test_loss2(unsigned int n, unsigned int m, unsigned int r, unsigned int c, unsigned int l, unsigned int a,
                       unsigned int N) // number of samples
{
    dcgp::kernel_set<gdual_d> basic_set({"sum", "diff", "mul", "div"});
    dcgp::expression<gdual_d> ex(n, m, r, c, l, a, basic_set(), 123);

    // creates N data points
    std::default_random_engine re;
    std::vector<std::vector<audi::gdual_d>> in;
    std::vector<std::vector<audi::gdual_d>> out;
    std::vector<audi::gdual_d> in_point(n);
    std::vector<audi::gdual_d> out_point(m);

    for (auto i = 0u; i < N; ++i) {
        // We only define the first node as a weight and we compute the derivative of the objfun wrt this.
        in_point[0] = audi::gdual_d(3, "w", 1);
        for (auto j = 1u; j < n; ++j) {
            in_point[j] = audi::gdual_d(std::uniform_real_distribution<double>(-1, 1)(re));
        }
        out_point = ex(in_point);
        for (auto &k : out_point) {
            k = audi::gdual_d(k.constant_cf());
        }
        in.push_back(in_point);
        out.push_back(out_point);
    }
    // computes the loss
    return ex.loss(in, out, "MSE", true);
}


BOOST_AUTO_TEST_CASE(construction)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    kernel_set<double> empty_set;
    std::vector<unsigned> arity{3, 2, 1};
    std::vector<unsigned> arity_wrong{3, 2, 1, 4};

    // Sanity checks
    BOOST_CHECK_THROW(expression<double>(0, 4, 2, 3, 4, arity, basic_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 0, 2, 3, 4, 2, basic_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 0, 3, 4, arity, basic_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 0, 4, 3, basic_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 0, arity, basic_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 4, 1, empty_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 4, arity_wrong, basic_set(), rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 4, {3, 0, 1}, basic_set(), rd()), std::invalid_argument);

    // Easy getters
    expression<double> ex(2, 4, 2, 3, 4, arity, basic_set(), rd());
    BOOST_CHECK_NO_THROW(ex.get());
    BOOST_CHECK_NO_THROW(ex({1.2, -0.34}));
    BOOST_CHECK_NO_THROW(ex({std::string("x"), std::string("y")}));
    BOOST_CHECK_NO_THROW(ex(std::vector<double>{1.2, -0.34}));
    BOOST_CHECK_NO_THROW(ex(std::vector<std::string>{"x", "y"}));
    BOOST_CHECK(ex.get_n() == 2u);
    BOOST_CHECK(ex.get_m() == 4u);
    BOOST_CHECK(ex.get_r() == 2u);
    BOOST_CHECK(ex.get_c() == 3u);
    BOOST_CHECK(ex.get_l() == 4u);
    BOOST_CHECK(ex.get_f().size() == 4u);
    BOOST_CHECK(ex.get_arity() == arity);

    // Checking validity of randomly produced chromosome
    BOOST_CHECK_NO_THROW(ex.is_valid(ex.get()));
}

BOOST_AUTO_TEST_CASE(compute)
{
    // Random seed
    std::random_device rd;

    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});

    /// Testing over Miller's test case from the PPSN 2014 tutorial
    expression<double> ex1(2, 4, 2, 3, 4, 2, basic_set(), rd());
    ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});

    CHECK_EQUAL_V(ex1({1., -1.}), std::vector<double>({0, -1, -1, 0}));
    CHECK_CLOSE_V(ex1({-.123, 2.345}), std::vector<double>({2.222, -0.288435, 0.676380075, 0}), 1e-8);

    ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 6, 5, 7, 3});
    CHECK_EQUAL_V(ex1({2., -3.}), std::vector<double>({6, -6, -18, 0}));
    CHECK_EQUAL_V(ex1({1., -1.}), std::vector<double>({2, -1, -1, 0}));
    CHECK_CLOSE_V(ex1({-.123, 2.345}), std::vector<double>({-4.69, -0.288435, 0.676380075, 0}), 1e-8);

    /// Testing over a single row program
    dcgp::expression<double> ex2(4, 1, 1, 10, 10, 2, basic_set(), rd());
    ex2.set({2, 3, 0, 0, 2, 2, 3, 0, 1, 1, 5,  4, 2, 6, 1, 0,
             7, 7, 3, 6, 7, 1, 7, 6, 2, 4, 10, 2, 3, 2, 10}); ///(x/y)/(2z-(t*x))
    CHECK_CLOSE_V(ex2({2., 3., 4., -2.}), std::vector<double>({0.055555555555555552}), 1e-8);
    CHECK_EQUAL_V(ex2({-1., 1., -1., 1.}), std::vector<double>({1}));
}

BOOST_AUTO_TEST_CASE(check_bounds)
{
    // Random seed
    std::random_device rd;

    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});

    // Testing an expression with arity 3, levels-back 2
    expression<double> ex(3, 1, 2, 3, 2, 3, basic_set(), rd());

    CHECK_EQUAL_V(ex.get_lb(), std::vector<unsigned int>(
                                   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 3, 3, 3, 5}));
    CHECK_EQUAL_V(ex.get_ub(), std::vector<unsigned int>(
                                   {3, 2, 2, 2, 3, 2, 2, 2, 3, 4, 4, 4, 3, 4, 4, 4, 3, 6, 6, 6, 3, 6, 6, 6, 8}));
}

BOOST_AUTO_TEST_CASE(get_gene_idx)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    expression<double> ex(3, 2, 3, 3, 2, {2, 1, 3}, basic_set(), rd());
    auto test = ex.get_gene_idx();
    BOOST_CHECK(test == (std::vector<unsigned>{0, 0, 0, 0, 3, 6, 9, 11, 13, 15, 19, 23}));
}

BOOST_AUTO_TEST_CASE(get_active_nodes_and_genes)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    expression<double> ex(3, 2, 3, 3, 4, {2, 1, 3}, basic_set(), rd());

    ex.set({0, 0, 1, 1, 1, 2, 2, 0, 2, 0, 3, 1, 2, 2, 4, 0, 6, 6, 7, 3, 6, 7, 8, 1, 7, 8, 2, 11, 11});

    auto a_nodes = ex.get_active_nodes();
    auto a_genes = ex.get_active_genes();
    auto node_x_idx = ex.get_gene_idx();
    for (auto i = 0u; i < a_nodes.size(); ++i) {
        auto node_id = a_nodes[i];
        // First the function gene
        if (node_id > ex.get_n()) {
            BOOST_CHECK(std::find(a_genes.begin(), a_genes.end(), node_x_idx[node_id]) != a_genes.end());
            // Then the connection genes
            for (auto j = 0u; j < ex.get_arity(node_id); ++j) {
                BOOST_CHECK(std::find(a_genes.begin(), a_genes.end(), node_x_idx[node_id] + j) != a_genes.end());
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(mutate)
{
    // Random seed
    std::random_device rd;
    // We define a d-CGP expression and get the reference chromosome
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    // Number of random trials
    unsigned int N = 100;
    // A d-CGP expression
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());

        // We test mutate(idx). Does the idx gene change? Do all other stay the
        // same?
        for (auto i = 0u; i < N; ++i) {
            std::vector<unsigned int> x = ex.get();
            unsigned idx = static_cast<unsigned>(rd() % x.size());
            ex.mutate(idx);
            for (auto j = 0u; j < x.size(); ++j) {
                if (j == idx) {
                    BOOST_CHECK(x[j] != ex.get()[j]);
                } else {
                    BOOST_CHECK_EQUAL(x[j], ex.get()[j]);
                }
            }
        }
    }
    // We test mutate_active. Was the mutated gene active?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i) {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_active();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned idx = static_cast<unsigned>(it.first - x.begin());
            auto ag = ex.get_active_genes();
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
        }
    }

    // We test mutate_active_fgene. Was the mutated gene active? Was it a function
    // gene?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i) {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_active_fgene();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned idx = static_cast<unsigned>(it.first - x.begin());
            auto ag = ex.get_active_genes();
            // was it an active gene?
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
            // was it a function gene?
            BOOST_CHECK(idx < x.size() - ex.get_m());
            BOOST_CHECK_EQUAL(idx % 3, 0);
        }
    }

    // We test mutate_active_cgene. Was the mutated gene active? Was it a
    // connection gene?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i) {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_active_cgene();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned idx = static_cast<unsigned>(it.first - x.begin());
            auto ag = ex.get_active_genes();
            // was it an active gene?
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
            // was it a connection gene?
            BOOST_CHECK(idx < x.size() - ex.get_m());
            BOOST_CHECK(((idx % 3) == 1 || (idx % 3) == 2) == true);
        }
    }

    // We test mutate_ogene. Was the mutated gene active? Was it an output gene?
    {
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), rd());
        for (auto i = 0u; i < N; ++i) {
            std::vector<unsigned int> x = ex.get();
            ex.mutate_ogene();
            auto it = std::mismatch(x.begin(), x.end(), ex.get().begin());
            unsigned idx = static_cast<unsigned>(it.first - x.begin());
            auto ag = ex.get_active_genes();
            // was it an active gene?
            BOOST_CHECK(std::find(ag.begin(), ag.end(), idx) != ag.end());
            // was it an output gene?
            BOOST_CHECK(idx >= x.size() - ex.get_m());
        }
    }
}

BOOST_AUTO_TEST_CASE(loss)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    expression<double> ex(2, 2, 2, 2, 3, 2, basic_set(), rd());
    // 2xy, 2x
    ex.set({0, 1, 1, 0, 0, 0, 2, 0, 2, 2, 0, 2, 4, 3});
    // MSE
    auto loss = ex.loss({1., 1.}, {2., 2.}, expression<double>::loss_type::MSE);
    BOOST_CHECK_EQUAL(loss, 0.);
    loss = ex.loss({1., 1.}, {0., 0.}, expression<double>::loss_type::MSE);
    BOOST_CHECK_EQUAL(loss, 4.);
    loss = ex.loss({1., 0.}, {0., 0.}, expression<double>::loss_type::MSE);
    BOOST_CHECK_EQUAL(loss, 2.);
    // CE
    loss = ex.loss({1., 1.}, {0.5, 0.5}, expression<double>::loss_type::CE);
    BOOST_CHECK_CLOSE(loss, 0.69314718055994529, 1e-12);
    loss = ex.loss({1., 1.}, {0.0, 1.0}, expression<double>::loss_type::CE);
    BOOST_CHECK_CLOSE(loss, 0.69314718055994529, 1e-12);
    loss = ex.loss({1., 0.}, {0., 1.}, expression<double>::loss_type::CE);
    BOOST_CHECK_CLOSE(loss, 0.12692801104297263, 1e-12);
    // On a batch (first sequential then parallel)
    auto loss_b = ex.loss({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, "MSE", false);
    BOOST_CHECK_EQUAL(loss_b, 1.);
    loss_b = ex.loss({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, "MSE", true);
    BOOST_CHECK_EQUAL(loss_b, 1.);
    // Identities
    // We test that a d-CGP expression computed on 20 points
    // has a zero quadratic error w.r.t. itself (its a perfect fit of itself)
    BOOST_CHECK_EQUAL(test_loss(3, 1, 1, 20, 21, 2, 20), 0);
    BOOST_CHECK_EQUAL(test_loss(2, 2, 3, 10, 11, 2, 20), 0);

    // We test that a d-CGP expression computed on 20 points
    // has a zero quadratic error w.r.t. itself, and that the
    // derivative of the quadratic error is zero w.r.t. one of the inputs (a weight)
    BOOST_CHECK_EQUAL(test_loss2(3, 1, 1, 20, 21, 2, 20), audi::gdual_d(0));
    BOOST_CHECK_EQUAL(test_loss2(2, 2, 3, 10, 11, 2, 20), audi::gdual_d(0));
}
