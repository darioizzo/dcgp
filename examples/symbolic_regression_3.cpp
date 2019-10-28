#include <boost/algorithm/string.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <symengine/expression.h>
#include <vector>

#include <dcgp/algorithms/mes4cgp.hpp>
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

    // We instantiate a symbolic regression problem with one only ephemeral constants.
    auto n_eph = 1u;
    symbolic_regression udp(X, Y, 1u, 15u, 16u, 2u, kernel_set<double>({"sum", "diff", "mul"})(), n_eph);

    // We init a population with four individuals.
    pagmo::population pop{udp, 4};

    // We instantiate the memetic solver setting 1000 maximum generation and one active mutation (minimum)
    dcgp::mes4cgp uda(1000u, 1u, 1e-8);
    pagmo::algorithm algo{uda};
    algo.set_verbosity(50u);

    // We solve
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
