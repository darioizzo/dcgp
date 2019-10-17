#include <boost/algorithm/string.hpp>
#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <symengine/expression.h>
#include <vector>

#include <dcgp/gym.hpp>

using namespace dcgp;
using namespace boost::algorithm;

// In this first tutorial we show how to find an exact formula for some input data that do not require
// any real valued constant.
//
// This is the easiest case for a symbolic regression task and thus makes it for a perfect entry tutorial.
//
// We use the problem P1 from the dcgp::gym, that is x0 - 2*x0**3 + x0**5.

int main()
{
    // We use the problem P1 in this example.
    std::vector<std::vector<double>> X, Y;
    gym::generate_koza_quintic(X, Y);

    // We instantiate a symbolic regression problem with no ephemeral constants
    symbolic_regression udp(X, Y, 1, 20, 21, 2, kernel_set<double>({"sum", "diff", "mul", "pdiv"})(), 0u);

    // We init a population with four individuals
    pagmo::population pop{udp, 4};

    // And we define an evolutionary startegy with 1000 generation and 2
    // active mutations (base)
    dcgp::es4cgp uda(10000, 2u, 1e-8);
    pagmo::algorithm algo{uda};
    algo.set_verbosity(100u);

    // We evolve the population
    pop = algo.evolve(pop);

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
