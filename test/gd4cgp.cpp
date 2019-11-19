#define BOOST_TEST_MODULE dcgp_es4cgp_test
#include <boost/test/included/unit_test.hpp>

#include <sstream>

#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/rosenbrock.hpp>

#include <dcgp/algorithms/gd4cgp.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <dcgp/s11n.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction_test)
{
    BOOST_CHECK_NO_THROW(gd4cgp(1u, 0.1, 0.001));
    BOOST_CHECK_THROW(gd4cgp(1u, 0.1, 0.2), std::invalid_argument);
    BOOST_CHECK_THROW(gd4cgp(1u, -0.1, 0.001), std::invalid_argument);
    BOOST_CHECK_THROW(gd4cgp(1u, 0.1, -0.001), std::invalid_argument);
    BOOST_CHECK_THROW(gd4cgp(1u, -0.1, -1.001), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(evolve_test)
{
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::vector<std::vector<double>> points, labels;
    gym::generate_koza_quintic(points, labels);
    // We test that the evolve fails on UDPs that are not suitable.
    gd4cgp uda(10u, 1u, 1e-4);
    { // wrong problem (not symbolic)
        pagmo::population pop{pagmo::rosenbrock(10), 4};
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // no eph constants
        pagmo::population pop(symbolic_regression(points, labels, 2, 2, 3, 2, basic_set(), 0u, 0u), 1u);
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // zero iters
        pagmo::population pop{symbolic_regression(points, labels, 2, 2, 3, 2, basic_set(), 2u, 0u), 4u};
        BOOST_CHECK(gd4cgp{0u}.evolve(pop).get_x()[0] == pop.get_x()[0]);
    }
}

BOOST_AUTO_TEST_CASE(trivial_methods_test)
{
    gd4cgp uda(10u, 1u, 1e-4);
    uda.set_verbosity(11u);
    BOOST_CHECK(uda.get_verbosity() == 11u);
    BOOST_CHECK(uda.get_name().find("gradient descent") != std::string::npos);
    BOOST_CHECK(uda.get_extra_info().find("Minimum learning rate") != std::string::npos);
    BOOST_CHECK_NO_THROW(uda.get_log());
}

BOOST_AUTO_TEST_CASE(pagmo_integration_test)
{
    // Here we test that the uda can be used to construct a pagmo algorithm. And we call one evolution
    // via pagmo to make sure the interface works.
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::vector<std::vector<double>> points, labels;
    gym::generate_koza_quintic(points, labels);
    gd4cgp uda(10u, 1u, 1e-4);
    BOOST_CHECK_NO_THROW(pagmo::algorithm{uda});
    pagmo::algorithm algo(uda);
    algo.set_verbosity(1u);
    pagmo::population pop(symbolic_regression(points, labels, 1, 16, 3, 2, basic_set(), 5u, 0u), 4u);
    pop = algo.evolve(pop);
}

BOOST_AUTO_TEST_CASE(s11n_test)
{
    gd4cgp uda(10u, 1u, 1e-4);

    const auto orig = uda.get_extra_info();

    std::stringstream ss;
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << uda;
    }
    uda = gd4cgp{10u, 2u, 1e-3};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> uda;
    }

    BOOST_CHECK(orig == uda.get_extra_info());
}
