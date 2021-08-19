#define BOOST_TEST_MODULE dcgp_mes4cgp_test
#include <boost/test/included/unit_test.hpp>

#include <sstream>

#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/rosenbrock.hpp>

#include <dcgp/gym.hpp>
#include <dcgp/algorithms/mes4cgp.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <dcgp/s11n.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction_test)
{
    BOOST_CHECK_NO_THROW(mes4cgp(0u, 1u, 0., 0u));
    BOOST_CHECK_THROW(mes4cgp(1u, 0u, 1e-4, 0u), std::invalid_argument);
    BOOST_CHECK_THROW(mes4cgp(1u, 1u, -1e-4, 0u), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(evolve_test)
{
    // We test that evolve fails on UDPs that are not suitable.
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    mes4cgp uda(10u, 1u, 1e-4, 0u);
    { // wrong problem (not symbolic)
        pagmo::population pop{pagmo::rosenbrock(10), 4};
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // small pop
        pagmo::population pop(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}), 1u);
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // multiobjective
        pagmo::population pop(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1, 15, 16, 2,
                                                  basic_set(), 2u, true),
                              10u);
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // zero gen
        pagmo::population pop{symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}), 4u};
        BOOST_CHECK(mes4cgp{0u}.evolve(pop).get_x()[0] == pop.get_x()[0]);
    }
    // Here we only test that evolution is deterministic if the seed is controlled
    pagmo::problem prob{symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}})};
    pagmo::population pop1{prob, 5u, 23u};
    pagmo::population pop2{prob, 5u, 23u};
    pagmo::population pop3{prob, 5u, 23u};

    mes4cgp uda1{10u, 1u, 1e-4, 23u};
    uda1.set_verbosity(1u);
    pop1 = uda1.evolve(pop1);
    BOOST_CHECK(uda1.get_log().size() > 0u);

    mes4cgp uda2{10u, 1u, 1e-4, 23u};
    uda2.set_verbosity(1u);
    pop2 = uda2.evolve(pop2);
    BOOST_CHECK(uda2.get_log() == uda1.get_log());

    uda2.set_seed(23u);
    pop3 = uda2.evolve(pop3);
    BOOST_CHECK(uda1.get_log() == uda2.get_log());
}

BOOST_AUTO_TEST_CASE(trivial_methods_test)
{
    mes4cgp uda{10u, 1u, 1e-4, 23u};
    uda.set_verbosity(11u);
    BOOST_CHECK(uda.get_verbosity() == 11u);
    uda.set_seed(5u);
    BOOST_CHECK(uda.get_seed() == 5u);
    BOOST_CHECK(uda.get_name().find("CGP") != std::string::npos);
    BOOST_CHECK(uda.get_extra_info().find("Verbosity") != std::string::npos);
    BOOST_CHECK_NO_THROW(uda.get_log());
}

BOOST_AUTO_TEST_CASE(s11n_test)
{
    mes4cgp uda{10u, 1u, 1e-4, 23u};

    const auto orig = uda.get_extra_info();

    std::stringstream ss;
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << uda;
    }
    uda = mes4cgp{10u, 1u, 1e-5, 23u};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> uda;
    }

    BOOST_CHECK(orig == uda.get_extra_info());
}

BOOST_AUTO_TEST_CASE(correct_fitness)
{
    mes4cgp uda{10u, 4u, 1e-4, 42u};

    kernel_set<double> basic_set({"sum", "diff", "mul", "pdiv"});
    std::vector<std::vector<double>> points, labels;
    gym::generate_koza_quintic(points, labels);
    symbolic_regression udp(points, labels, 1u, 10u, 11u, 2u, basic_set(), 3u, false, 0u, "MSE", 42u);
    
    pagmo::population pop1(udp, 4u, 16u);
    
    pop1 = uda.evolve(pop1);
    // evolve population once 

    pagmo::population pop2(udp, 0u);

    pop2.push_back(pop1.champion_x());
    BOOST_CHECK(pop2.champion_f() == pop1.champion_f());
}
