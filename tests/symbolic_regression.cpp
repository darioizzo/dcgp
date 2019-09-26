#define BOOST_TEST_MODULE dcgp_symbolic_regression_test
#include <boost/test/unit_test.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/sga.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction_test)
{
    // Its default-onstructable
    BOOST_CHECK_NO_THROW(symbolic_regression{});
    // Sanity checks tests (inconsistent points / labels)
    BOOST_CHECK_THROW(symbolic_regression({}, {}), std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}, {0.02 / 0.32}}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32, 0.3}}, {{3. / 2.}, {0.02 / 0.32}}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2., 2.2}, {0.02 / 0.32}}),
                      std::invalid_argument);
    // Sanity checks tests (inconsistent cgp parameters)
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 0u, 1u, 1u, 2u,
                                          kernel_set<double>({"sum"})(), 0u),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1u, 0u, 1u, 2u,
                                          kernel_set<double>({"sum"})(), 0u),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1u, 1u, 0u, 2u,
                                          kernel_set<double>({"sum"})(), 0u),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1u, 1u, 1u, 1u,
                                          kernel_set<double>({"sum"})(), 0u),
                      std::invalid_argument);
    BOOST_CHECK_THROW(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1u, 1u, 1u, 2u,
                                          kernel_set<double>(std::vector<std::string>())(), 0u),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(fitness_test)
{
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    // 2xy, 2x
    pagmo::vector_double test_x = {0, 1, 1, 0, 0, 0, 2, 0, 2, 2, 0, 2, 4, 3};
    // On a single point/label.
    {
        symbolic_regression udp({{1., 1.}}, {{2., 2.}}, 2, 2, 3, 2, basic_set(), 0u);
        BOOST_CHECK_EQUAL(udp.fitness(test_x)[0], 0.);
    }
    {
        symbolic_regression udp({{1., 1.}}, {{0., 0.}}, 2, 2, 3, 2, basic_set(), 0u);
        BOOST_CHECK_EQUAL(udp.fitness(test_x)[0], 4.);
    }
    {
        symbolic_regression udp({{1., 0.}}, {{0., 0.}}, 2, 2, 3, 2, basic_set(), 0u);
        BOOST_CHECK_EQUAL(udp.fitness(test_x)[0], 2.);
    }
    // On a batch (first sequential then parallel)
    {
        symbolic_regression udp({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, 2, 2, 3, 2, basic_set(), 0u);
        BOOST_CHECK_EQUAL(udp.fitness(test_x)[0], 1.);
    }
    {
        symbolic_regression udp({{1., 1.}, {1., 0.}}, {{2., 2.}, {0., 0.}}, 2, 2, 3, 2, basic_set(), 1u);
        BOOST_CHECK_EQUAL(udp.fitness(test_x)[0], 1.);
    }
}

BOOST_AUTO_TEST_CASE(get_bounds_test)
{
    symbolic_regression udp({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}});
    auto cgp = udp.get_cgp();
    auto lbu = cgp.get_lb();
    std::vector<double> lb(lbu.size());
    std::transform(lbu.begin(), lbu.end(), lb.begin(), [](unsigned a) { return boost::numeric_cast<double>(a); });
    BOOST_CHECK(udp.get_bounds().first == lb);
}

BOOST_AUTO_TEST_CASE(trivial_methods_test)
{
    symbolic_regression udp({{1., 1.}}, {{2., 2.}}, 2, 2, 3, 2);
    BOOST_CHECK_EQUAL(udp.get_bounds().first.size(), udp.get_nix());
    BOOST_CHECK(udp.get_name().find("CGP") != std::string::npos);
    BOOST_CHECK(udp.get_extra_info().find("Input dimension") != std::string::npos);
    pagmo::vector_double test_x = {0, 1, 1, 0, 0, 0, 2, 0, 2, 2, 0, 2, 4, 3};
    BOOST_CHECK(udp.pretty(test_x).find("[(x0*(x1+x1)), (x0+x0)]") != std::string::npos);
    BOOST_CHECK_NO_THROW(udp.get_cgp());
}