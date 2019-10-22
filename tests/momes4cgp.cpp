#define BOOST_TEST_MODULE dcgp_momes4cgp_test
#include <boost/test/unit_test.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/rosenbrock.hpp>

#include <dcgp/algorithms/momes4cgp.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <dcgp/gym.hpp>


using namespace dcgp;

BOOST_AUTO_TEST_CASE(evolve_test)
{
    // The Data
    std::vector<std::vector<double>> points, labels;
    gym::generate_vladi4(points, labels);
    // The problem (multi-objective)
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    symbolic_regression udp{points, labels, 1, 15, 16, 2, basic_set(), 2, true, 0u};
    // The algorithm 
    momes4cgp uda{100, 4};
    uda.set_verbosity(1);
    // The initial population
    pagmo::population pop{udp, 20};
    //The Evolution
    pop = uda.evolve(pop);


}