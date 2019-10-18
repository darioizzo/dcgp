#define BOOST_TEST_MODULE dcgp_es4cgp_test
#include <boost/test/unit_test.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/rosenbrock.hpp>

#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction_test)
{
    BOOST_CHECK_NO_THROW(es4cgp(0u, 2u, 0., true, 0u));
    BOOST_CHECK_THROW(es4cgp(1u, 0u, 1e-4, true, 0u), std::invalid_argument);
    BOOST_CHECK_THROW(es4cgp(1u, 2u, -1e-4, true, 0u), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(evolve_test)
{
    // We test that evolve fails on UDPs that are not suitable.
    es4cgp uda(10u, 2u, 1e-4, true, 0u);
    { // wrong problem (not symbolic)
        pagmo::population pop{pagmo::rosenbrock(10), 4};
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // small pop
        pagmo::population pop(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}), 1u);
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // zero gen
        pagmo::population pop{symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}), 4u};
        BOOST_CHECK(es4cgp{0u}.evolve(pop).get_x()[0] == pop.get_x()[0]);
    }
    // Here we only test that evolution is deterministic if the
    // seed is controlled for all variants
    pagmo::problem prob{symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}})};
    pagmo::population pop1{prob, 5u, 23u};
    pagmo::population pop2{prob, 5u, 23u};
    pagmo::population pop3{prob, 5u, 23u};

    es4cgp uda1{10u, 2u, 1e-4, true, 23u};
    uda1.set_verbosity(1u);
    pop1 = uda1.evolve(pop1);
    BOOST_CHECK(uda1.get_log().size() > 0u);

    es4cgp uda2{10u, 2u, 1e-4, true, 23u};
    uda2.set_verbosity(1u);
    pop2 = uda2.evolve(pop2);
    BOOST_CHECK(uda2.get_log() == uda1.get_log());

    uda2.set_seed(23u);
    pop3 = uda2.evolve(pop3);
    BOOST_CHECK(uda1.get_log() == uda2.get_log());
}

BOOST_AUTO_TEST_CASE(trivial_methods_test)
{
    es4cgp uda{10u, 2u, 1e-4, true, 23u};
    uda.set_verbosity(11u);
    BOOST_CHECK(uda.get_verbosity() == 11u);
    uda.set_seed(5u);
    BOOST_CHECK(uda.get_seed() == 5u);
    BOOST_CHECK(uda.get_name().find("CGP") != std::string::npos);
    BOOST_CHECK(uda.get_extra_info().find("Verbosity") != std::string::npos);
    BOOST_CHECK_NO_THROW(uda.get_log());
}