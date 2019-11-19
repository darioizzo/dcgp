#define BOOST_TEST_MODULE dcgp_momes4cgp_test
#include <boost/test/included/unit_test.hpp>

#include <sstream>

#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/rosenbrock.hpp>

#include <dcgp/algorithms/momes4cgp.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <dcgp/s11n.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction_test)
{
    BOOST_CHECK_NO_THROW(momes4cgp(0u, 1u, 0u));
    BOOST_CHECK_THROW(momes4cgp(1u, 0u, 0u), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(evolve_test)
{
    // We test that evolve fails on UDPs that are not suitable.
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    momes4cgp uda(10u, 1u, 0u);
    { // wrong problem (not symbolic)
        pagmo::population pop{pagmo::rosenbrock(10), 4};
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // small pop
        pagmo::population pop(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}), 1u);
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // singleobjective
        pagmo::population pop(symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1, 15, 16, 2,
                                                  basic_set(), 2u, false),
                              10u);
        BOOST_CHECK_THROW(uda.evolve(pop), std::invalid_argument);
    }
    { // zero gen
        pagmo::population pop{symbolic_regression({{1., 2.}, {0.3, -0.32}}, {{3. / 2.}, {0.02 / 0.32}}, 1, 20, 21, 2,
                                                  basic_set(), 1u, true, 0u),
                              4u};
        BOOST_CHECK(momes4cgp{0u}.evolve(pop).get_x()[0] == pop.get_x()[0]);
    }
    // Here we only test that evolution is deterministic if the seed is controlled
    std::vector<std::vector<double>> points, labels;
    gym::generate_koza_quintic(points, labels);
    pagmo::problem prob{symbolic_regression(points, labels, 1, 20, 21, 2, basic_set(), 1u, true, 0u)};

    pagmo::population pop1{prob, 5u, 23u};
    pagmo::population pop2{prob, 5u, 23u};
    pagmo::population pop3{prob, 5u, 23u};

    momes4cgp uda1{10u, 1u, 23u};
    uda1.set_verbosity(1u);
    pop1 = uda1.evolve(pop1);
    BOOST_CHECK(uda1.get_log().size() > 0u);

    momes4cgp uda2{10u, 1u, 23u};
    uda2.set_verbosity(1u);
    pop2 = uda2.evolve(pop2);
    BOOST_CHECK(uda2.get_log() == uda1.get_log());

    uda2.set_seed(23u);
    pop3 = uda2.evolve(pop3);
    BOOST_CHECK(uda1.get_log() == uda2.get_log());
}

BOOST_AUTO_TEST_CASE(trivial_methods_test)
{
    momes4cgp uda{10u, 1u, 23u};
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
    momes4cgp uda{10u, 1u, 23u};

    const auto orig = uda.get_extra_info();

    std::stringstream ss;
    {
        boost::archive::binary_oarchive oarchive(ss);
        oarchive << uda;
    }
    uda = momes4cgp{10u, 2u, 23u};
    {
        boost::archive::binary_iarchive iarchive(ss);
        iarchive >> uda;
    }

    BOOST_CHECK(orig == uda.get_extra_info());
}
