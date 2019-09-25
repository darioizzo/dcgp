#include <dcgp/algorithms/es4cgp.hpp>
#include <dcgp/problems/symbolic_regression.hpp>
#include <pagmo/algorithm.hpp>
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
    pagmo::population pop{prob, 4};
    dcgp::es4cgp uda(10000);
    pagmo::algorithm algo{uda};
    pagmo::print(algo);
    algo.set_verbosity(100u);
    pop = algo.evolve(pop);
    pagmo::print(pop);
    return false;
}
