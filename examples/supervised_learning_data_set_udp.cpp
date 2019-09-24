#include <dcgp/problems/symbolic_regression.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>
#include <pagmo/algorithms/sga.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <vector>

// reading a text file
#include "detail/read_data.hpp"
using namespace dcgp;

int main()
{
    // We read the data from file
    std::vector<std::vector<double>> X, Y;
    read_data(X, Y, "../../examples/data/symbolic.data");
    symbolic_regression udp(X, Y);
    pagmo::problem prob{udp};
    pagmo::print(prob);
    pagmo::population pop{prob, 10};
    // A Genetic algorithm without crossover and with only mutation
    pagmo::sga uda(500, 0., 1., 0.1);
    pagmo::algorithm algo{uda};
    pagmo::print(algo);
    algo.set_verbosity(1u);
    pop = algo.evolve(pop);
    // We print out the final expression
    audi::print("Final expression: ",  udp.pretty(pop.champion_x()), "\n");
    audi::print("Final value: ", pop.get_f()[pop.best_idx()], "\n");
    return false;
}
