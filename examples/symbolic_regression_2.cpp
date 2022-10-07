#include <boost/algorithm/string.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <symengine/expression.h>
// patch to avoid flint defining access come _access and messing with boost
#if defined(access)
#undef access
#endif
#include <vector>

#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;
using namespace boost::algorithm;

int main()
{
    // We load the data (using problem P1 from the gym)
    std::vector<std::vector<double>> X, Y;
    gym::generate_P1(X, Y);

    // We instantiate a symbolic regression problem with one ephemeral constants, and multiobjectove.
    auto n_eph = 1u;
    symbolic_regression udp(X, Y, 1u, 20u, 21u, 2u, kernel_set<double>({"sum", "diff", "mul", "pdiv"})(), n_eph, false,
                            0u, "MSE");

    // We init a population with four individuals
    pagmo::population pop{udp, 4u};

    // We use an ES startegy also to learn constants
    dcgp::es4cgp uda1(100000u, 4u, 1e-8, true);
    pagmo::algorithm algo{uda1};
    algo.set_verbosity(500u);

    auto t1 = std::chrono::high_resolution_clock::now();
    pop = algo.evolve(pop);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << duration;

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
