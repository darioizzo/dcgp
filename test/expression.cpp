#define BOOST_TEST_MODULE dcgp_compute_test
#include <boost/test/included/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <audi/gdual.hpp>

#include <pagmo/io.hpp>

#include <dcgp/expression.hpp>
#include <dcgp/function.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/s11n.hpp>
#include <dcgp/wrapped_functions_s11n_implement.hpp>

#include "helpers.hpp"

using namespace dcgp;
using namespace audi;

double test_loss(unsigned int n, unsigned int m, unsigned int r, unsigned int c, unsigned int l, unsigned int a,
                 unsigned int N) // number of samples
{
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    dcgp::expression<double> ex(n, m, r, c, l, a, basic_set(), 0u, 123u);

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

gdual_d test_loss2(unsigned int n, unsigned int m, unsigned int r, unsigned int c, unsigned int l, unsigned int a,
                   unsigned int N) // number of samples
{
    dcgp::kernel_set<gdual_d> basic_set({"sum", "diff", "mul", "div"});
    dcgp::expression<gdual_d> ex(n, m, r, c, l, a, basic_set(), 0u, 123u);

    // creates N data points
    std::default_random_engine re;
    std::vector<std::vector<gdual_d>> in;
    std::vector<std::vector<gdual_d>> out;
    std::vector<gdual_d> in_point(n);
    std::vector<gdual_d> out_point(m);

    for (auto i = 0u; i < N; ++i) {
        // We only define the first node as a weight and we compute the derivative of the objfun wrt this.
        in_point[0] = gdual_d(3, "w", 1);
        for (auto j = 1u; j < n; ++j) {
            in_point[j] = gdual_d(std::uniform_real_distribution<double>(-1, 1)(re));
        }
        out_point = ex(in_point);
        for (auto &k : out_point) {
            k = gdual_d(k.constant_cf());
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
    BOOST_CHECK_THROW(expression<double>(0, 4, 2, 3, 4, arity, basic_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 0, 2, 3, 4, 2, basic_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 0, 3, 4, arity, basic_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 0, 4, 3, basic_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 0, arity, basic_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 4, 1, empty_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 4, arity_wrong, basic_set(), 0u, rd()), std::invalid_argument);
    BOOST_CHECK_THROW(expression<double>(2, 4, 2, 3, 4, {3, 0, 1}, basic_set(), 0u, rd()), std::invalid_argument);

    // Easy getters
    expression<double> ex(2, 4, 2, 3, 4, arity, basic_set(), 0u, rd());
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
}

BOOST_AUTO_TEST_CASE(compute)
{
    // Random seed
    std::random_device rd;

    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});

    /// Testing over Miller's test case from the PPSN 2014 tutorial
    expression<double> ex1(2, 4, 2, 3, 4, 2, basic_set(), 0u, rd());
    ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});

    CHECK_EQUAL_V(ex1({1., -1.}), std::vector<double>({0, -1, -1, 0}));
    CHECK_CLOSE_V(ex1({-.123, 2.345}), std::vector<double>({2.222, -0.288435, 0.676380075, 0}), 1e-8);

    ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 6, 5, 7, 3});
    CHECK_EQUAL_V(ex1({2., -3.}), std::vector<double>({6, -6, -18, 0}));
    CHECK_EQUAL_V(ex1({1., -1.}), std::vector<double>({2, -1, -1, 0}));
    CHECK_CLOSE_V(ex1({-.123, 2.345}), std::vector<double>({-4.69, -0.288435, 0.676380075, 0}), 1e-8);

    /// Testing over a single row program
    dcgp::expression<double> ex2(4, 1, 1, 10, 10, 2, basic_set(), 0u, rd());
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
    expression<double> ex(3, 1, 2, 3, 2, 3, basic_set(), 0u, rd());

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
    expression<double> ex(3, 2, 3, 3, 2, {2, 1, 3}, basic_set(), 0u, rd());
    auto test = ex.get_gene_idx();
    BOOST_CHECK(test == (std::vector<unsigned>{0, 0, 0, 0, 3, 6, 9, 11, 13, 15, 19, 23}));
}

BOOST_AUTO_TEST_CASE(get_active_nodes_and_genes)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    expression<double> ex(3, 2, 3, 3, 4, {2, 1, 3}, basic_set(), 0u, rd());

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
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), 0u, rd());

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
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), 0u, rd());
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
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), 0u, rd());
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
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), 0u, rd());
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
        expression<double> ex(3, 3, 2, 20, 21, 2, basic_set(), 0u, rd());
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
    expression<double> ex(2, 2, 2, 2, 3, 2, basic_set(), 0u, rd());
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
    auto loss_b = ex.loss({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, "MSE", 0u);
    BOOST_CHECK_EQUAL(loss_b, 1.);
    loss_b = ex.loss({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, "MSE", 1u);
    BOOST_CHECK_EQUAL(loss_b, 1.);
    // Identities
    // We test that a d-CGP expression computed on 20 points
    // has a zero quadratic error w.r.t. itself (its a perfect fit of itself)
    BOOST_CHECK_EQUAL(test_loss(3, 1, 1, 20, 21, 2, 20), 0);
    BOOST_CHECK_EQUAL(test_loss(2, 2, 3, 10, 11, 2, 20), 0);

    // We test that a d-CGP expression computed on 20 points
    // has a zero quadratic error w.r.t. itself, and that the
    // derivative of the quadratic error is zero w.r.t. one of the inputs (a weight)
    BOOST_CHECK_EQUAL(test_loss2(3, 1, 1, 20, 21, 2, 20), gdual_d(0));
    BOOST_CHECK_EQUAL(test_loss2(2, 2, 3, 10, 11, 2, 20), gdual_d(0));

    std::mt19937 mersenne_engine{rd()}; // Generates random integers
    std::uniform_real_distribution<double> dist{-1., 1.};

    // We test that the parallel and the sequantial algorithms both return the same result
    for (auto i = 0u; i < 100; ++i) {
        auto in = std::vector<std::vector<double>>(100, {0., 0.});
        auto out = std::vector<std::vector<double>>(100, {0., 0.});
        std::generate(in.begin(), in.end(), [&mersenne_engine, &dist]() {
            return std::vector<double>{dist(mersenne_engine), dist(mersenne_engine)};
        });
        std::generate(out.begin(), out.end(), [&mersenne_engine, &dist]() {
            return std::vector<double>{dist(mersenne_engine), dist(mersenne_engine)};
        });
        BOOST_CHECK_CLOSE(ex.loss(in, out, "MSE", true), ex.loss(in, out, "MSE", false), 1e-8);
        BOOST_CHECK_CLOSE(ex.loss(in, out, "CE", true), ex.loss(in, out, "CE", false), 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(ephemeral_constants_test)
{
    std::vector<unsigned> test_x
        = {3, 1, 2, 2, 1, 4, 3, 1, 1, 3, 6, 5, 0, 0, 3, 1, 8, 0, 1, 2, 5, 0, 0, 8, 1, 3, 8, 0, 11, 11, 4, 12};
    {
        kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
        expression<double> ex(2, 2, 1, 10, 11, 2, basic_set(), 2u, 123u);
        // [(y/c1), (c2-(x+c2))]

        ex.set(test_x);
        BOOST_CHECK(ex({"x", "y"})[0].compare("(y/c1)") == 0);
        BOOST_CHECK(ex({"x", "y"})[1].compare("(c2-(x+c2))") == 0);

        // We check the setter/getter for ephemeral constants names
        BOOST_CHECK_THROW(ex.set_eph_symb({"C1", "C2", "C3"}), std::invalid_argument);
        ex.set_eph_symb({"C1", "C2"});
        BOOST_CHECK(ex({"x", "y"})[0].compare("(y/C1)") == 0);
        BOOST_CHECK(ex({"x", "y"})[1].compare("(C2-(x+C2))") == 0);
        BOOST_CHECK((ex.get_eph_symb() == std::vector<std::string>{"C1", "C2"}));

        // We check the setter/getter for ephemeral constants values
        BOOST_CHECK_THROW(ex.set_eph_val({1., 2., 3.}), std::invalid_argument);
        ex.set_eph_val({1., 2.});
        BOOST_CHECK((ex({3., 4.}) == std::vector<double>{4., -3.}));
        ex.set_eph_val({2., 0.});
        BOOST_CHECK((ex({3., 4.}) == std::vector<double>{2., -3.}));
        BOOST_CHECK((ex.get_eph_val() == std::vector<double>{2., 0.}));
    }
    // More on Setters and getters
    {
        kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
        expression<double> ex(2, 2, 1, 10, 11, 2, basic_set(), 2u, 123u);
        ex.set_eph_val({3.1, -1.2});
        BOOST_CHECK((ex.get_eph_val() == std::vector<double>{3.1, -1.2}));
        BOOST_CHECK((ex.get_eph_symb() == std::vector<std::string>{"c1", "c2"}));
        ex.set_eph_symb({"x", "y"});
        BOOST_CHECK((ex.get_eph_symb() == std::vector<std::string>{"x", "y"}));
    }
    // Lastly we check (a bit randomly) that this works with gduals and a loss (the basis for grad descent)
    {
        kernel_set<gdual_d> basic_set({"sum", "diff", "mul", "div"});
        expression<gdual_d> ex(2, 2, 1, 10, 11, 2, basic_set(), 2u, 123u);
        ex.set(test_x);
        ex.set_eph_val({gdual_d(1., "c1", 3), gdual_d(2., "c2", 3)});
        gdual_d c1(0., "c1", 3);
        BOOST_CHECK_EQUAL(4. - 4. * c1 + 4 * c1 * c1 - 4 * c1 * c1 * c1, ex({gdual_d(3.), gdual_d(4.)})[0]);
        BOOST_CHECK_EQUAL(gdual_d(-3), ex({gdual_d(3.), gdual_d(4.)})[1]);
        std::vector<std::vector<gdual_d>> in = {{gdual_d(-1.), gdual_d(1.)}, {gdual_d(-0.5), gdual_d(0.01)}};
        std::vector<std::vector<gdual_d>> out = {{gdual_d(0.), gdual_d(0.)}, {gdual_d(-2.), gdual_d(3.)}};
        BOOST_CHECK_CLOSE(ex.loss(in, out, "MSE", true).get_derivative({{"dc1", 1u}}), -0.51005, 1e-3);
    }
}

template <typename T>
std::vector<T> my_pc(const std::vector<T> &x, dcgp::function<std::vector<T>(const std::vector<T>&)>g_f)
{
    std::vector<T> retval(4, T(0));
    auto g = g_f(x);
    retval[0] = g[0] * x[0];
    retval[1] = g[1] * x[0] * x[1];
    retval[2] = g[2] * x[1];
    retval[3] = g[3] + x[0] * x[1];
    return retval;
}

template <typename T>
std::vector<T> my_pc2(const std::vector<T> &x, dcgp::function<std::vector<T>(const std::vector<T>&)>g_f)
{
    std::vector<T> retval(1, T(0));
    auto g = g_f(x);
    retval[0] = g[0] * x[0] * x[1] * x[2] * x[3];
    return retval;
}

BOOST_AUTO_TEST_CASE(phenotype_correction)
{
    // Random seed
    std::random_device rd;
    using pc_fun_type = expression<double>::pc_fun_type;
    // testing on double
    {
        kernel_set<double> basic_set({"sum", "diff", "mul", "div"});

        /// Testing over Miller's test case from the PPSN 2014 tutorial
        expression<double> ex1(2, 4, 2, 3, 4, 2, basic_set(), 0u, rd());
        // Test setter and values
        ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});
        ex1.set_phenotype_correction(pc_fun_type(my_pc<double>));
        CHECK_EQUAL_V(ex1({1., -1.}), std::vector<double>({0 * 1., -1 * 1. * -1., -1. * -1., 0 + -1.}));
        CHECK_CLOSE_V(
            ex1({-.123, 2.345}),
            std::vector<double>({2.222 * -.123, -0.288435 * -.123 * 2.345, 0.676380075 * 2.345, 0 - 2.345 * 0.123}),
            1e-8);
        // Test the unsetter
        ex1.unset_phenotype_correction();
        CHECK_EQUAL_V(ex1({1., -1.}), std::vector<double>({0, -1, -1, 0}));
        CHECK_CLOSE_V(ex1({-.123, 2.345}), std::vector<double>({2.222, -0.288435, 0.676380075, 0}), 1e-8);

        /// Testing over a single row program
        dcgp::expression<double> ex2(4, 1, 1, 10, 10, 2, basic_set(), 0u, rd());
        ex2.set({2, 3, 0, 0, 2, 2, 3, 0, 1, 1, 5,  4, 2, 6, 1, 0,
                 7, 7, 3, 6, 7, 1, 7, 6, 2, 4, 10, 2, 3, 2, 10}); ///(x/y)/(2z-(t*x))
        ex2.set_phenotype_correction(pc_fun_type(my_pc2<double>));
        CHECK_CLOSE_V(ex2({2., 3., 4., -2.}), std::vector<double>({0.055555555555555552 * 2. * 3. * 4. * -2.}), 1e-8);
        CHECK_EQUAL_V(ex2({-1., 1., -1., 1.}), std::vector<double>({1. * -1. * 1. * -1. * 1.}));
        ex2.unset_phenotype_correction();
        CHECK_CLOSE_V(ex2({2., 3., 4., -2.}), std::vector<double>({0.055555555555555552}), 1e-8);
        CHECK_EQUAL_V(ex2({-1., 1., -1., 1.}), std::vector<double>({1}));
    }
    // testing on gduals
    {
        kernel_set<gdual_d> basic_set({"sum", "diff", "mul", "div"});

        /// Testing over Miller's test case from the PPSN 2014 tutorial
        expression<gdual_d> ex1(2, 4, 2, 3, 4, 2, basic_set(), 0u, rd());
        // Test setter and values
        ex1.set({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});
        ex1.set_phenotype_correction(pc_fun_type(my_pc<gdual_d>));
        gdual_d x0(1.234, "x0", 1);
        gdual_d x1(-1.234, "x1", 1);
        // f(x, g(x))
        auto f_g = ex1({x0, x1});
        ex1.unset_phenotype_correction();
        // g(x)
        auto g_g = ex1({x0, x1});
        auto f = f_g[0].constant_cf();
        auto g = g_g[0].constant_cf();
        BOOST_CHECK_EQUAL(f, g * 1.234);
        auto df = f_g[0].get_derivative({{"dx0", 1u}});
        auto dg = g_g[0].get_derivative({{"dx0", 1u}});
        BOOST_CHECK_EQUAL(df, dg * 1.234 + g);
    }
}

struct my_pc3 {
    /// Call operator
    std::vector<double> operator()(const std::vector<double> &x, dcgp::function<std::vector<double>(const std::vector<double>&)>g_f) const
    {
        std::vector<double> retval = g_f(x);
        retval[0] = retval[0]*2;
        retval[1] = retval[1]*10;
        return retval;
    }
    DCGP_S11N_EMPTY_SERIALIZE_MEMFN()
};
DCGP_S11N_FUNCTION_EXPORT_KEY(my_pc3_double, my_pc3, std::vector<double>, const std::vector<double> &, dcgp::function<std::vector<double>(const std::vector<double>&)>)
DCGP_S11N_FUNCTION_IMPLEMENT(my_pc3_double, my_pc3, std::vector<double>, const std::vector<double> &, dcgp::function<std::vector<double>(const std::vector<double>&)>)



BOOST_AUTO_TEST_CASE(s11n_test)
{
    // Random seed
    std::random_device rd;
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    expression<double> ex(2, 2, 2, 2, 3, 2, basic_set(), 0u, rd());
    ex.set_phenotype_correction(my_pc3());
    auto before_num = ex({1.2, 3.3});
    auto before_string = boost::lexical_cast<std::string>(ex);

    std::stringstream ss;
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << ex;
    }
    ex = expression<double>(2, 2, 2, 2, 3, 2, basic_set(), 0u, rd());
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> ex;
    }
    auto after_num = ex({1.2, 3.3});

    BOOST_CHECK(before_num == after_num);
    BOOST_CHECK(before_string == boost::lexical_cast<std::string>(ex));

}
