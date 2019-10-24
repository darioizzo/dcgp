#include <boost/algorithm/string.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <symengine/expression.h>
#include <vector>

#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;
using namespace boost::algorithm;

// In this second tutorial we show how to find a formula for some input data whn we also
// want to learn some constants.
//
// This case is more generic and interesting than the one treated in the previous tutorial,
// and can be considered, already, of use in industrial settings.
//
// NOTE: constants can, in general, be learned via two main techniques
// 1) evolutionary (common practice in GP)
// 2) memetic (this is original with dCGP)
//
// In this tutorial we follow the evolutionary approach 1). In the next tutorial we will follow a memetic approach 2).
//
// We use the problem P1 from the dcgp::gym, that is x**5 - pi*x**3 + x

int main()
{
    // We load the data (using problem P1 from the gym)
    std::vector<std::vector<double>> X, Y;
    gym::generate_P1(X, Y);

    // We instantiate a symbolic regression problem with one ephemeral constants
    auto n_eph = 1u;
    symbolic_regression udp(X, Y, 1u, 15u, 16u, 2u, kernel_set<double>({"sum", "diff", "mul"})(), n_eph);

    // We init a population with four individuals
    pagmo::population pop{udp, 4};

    // We use an ES startegy also to learn constants
    dcgp::es4cgp uda1(10000, 1u, 1e-8, true);
    pagmo::algorithm algo_gobal{uda1};
    algo_gobal.set_verbosity(500u);
    pop = algo_gobal.evolve(pop);

    // We print on screen the best found
    auto idx = pop.best_idx();
    auto prettier = udp.prettier(pop.get_x()[idx]);
    trim_left_if(prettier, is_any_of("["));
    trim_right_if(prettier, is_any_of("]"));

    pagmo::print("\nBest fitness: ", pop.get_f()[idx], "\n");
    pagmo::print("Chromosome: ", pop.get_x()[idx], "\n");
    pagmo::print("Pretty Formula: ", udp.pretty(pop.get_x()[idx]), "\n");
    pagmo::print("Prettier Formula: ", prettier, "\n");
    pagmo::print("Expanded Formula: ", SymEngine::expand(SymEngine::Expression(prettier)), "\n");
    return false;
}
