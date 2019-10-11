#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <vector>

#include <dcgp/gym.hpp>

// reading a text file
#include "detail/read_data.hpp"

using namespace dcgp;

int main()
{
    // We read the data from file
    std::vector<std::vector<double>> X, Y;
    gym::generate_P1(X, Y);
    // read_data(X, Y, "../../examples/data/symbolic.data");
    symbolic_regression udp(X, Y, 1, 20, 21, 2, kernel_set<double>({"sum", "diff", "mul", "pdiv"})(), 1u);
    pagmo::problem prob{udp};
    pagmo::population pop{prob, 4};
    dcgp::es4cgp uda(10000, 2u, 1e-8);
    pagmo::algorithm algo{uda};
    algo.set_verbosity(100u);
    pop = algo.evolve(pop);
    auto idx = pop.best_idx();
    pagmo::print(pop.get_f()[idx], "\n");
    pagmo::print(udp.pretty(pop.get_x()[idx]), "\n");
    return false;
}
