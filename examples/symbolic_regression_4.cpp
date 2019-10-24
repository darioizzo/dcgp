#include <boost/algorithm/string.hpp>
#include <pagmo/algorithm.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/utils/multi_objective.hpp>
#include <symengine/expression.h>
#include <vector>

#include <dcgp/algorithms/momes4cgp.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;
using namespace boost::algorithm;

// In this fourth tutorial we solve the multiobjective symbolic regression problem
// where the loss is considered next to the formula complexity to determine how good 
// a certain expression is.
//
// We use the problem vladi6 from the dcgp::gym, that is 6*cos(x0*sin(x1))

int main()
{
    // We load the data (using problem vladi6 from the gym)
    std::vector<std::vector<double>> X, Y;
    gym::generate_vladi6(X, Y);

    // We instantiate a symbolic regression problem with one only ephemeral constants.
    // Note that here we also set a batch parallelism to 5 so that 5 batches of 2 points
    // will be created and handled in parallel.
    auto n_eph = 1u;
    symbolic_regression udp(X, Y, 1u, 100, 101, 2u, kernel_set<double>({"sum", "diff", "mul", "sin", "cos"})(), n_eph, true, 5u);

    // We init a large population (100) of individuals.
    pagmo::population pop{udp, 100};

    // We here use the Multi-Objective Memetic Evolutionary Strategy: an original algorithm provided in this dCGP project.
    // We instantiate it with 100 generations and 4 active mutations maximum per individual.
    dcgp::momes4cgp uda{500u, 4u};
    pagmo::algorithm algo{uda};
    algo.set_verbosity(1u);

    // We solve the problem
    pop = algo.evolve(pop);

    // Finally we print on screen the non dominated front.
    pagmo::print("\n");
    auto ndf = pagmo::non_dominated_front_2d(pop.get_f());
    for (decltype(ndf.size()) i = 0u; i < ndf.size(); ++i) {
        auto idx = ndf[i];
        auto prettier = udp.prettier(pop.get_x()[idx]);
        trim_left_if(prettier, is_any_of("["));
        trim_right_if(prettier, is_any_of("]"));
        pagmo::print("Loss: ", std::setw(10), std::left, pop.get_f()[idx][0], std::setw(15), "Complexity: ", std::left, std::setw(5), pop.get_f()[idx][1], std::setw(10), "Formula: ", prettier,"\n");
    }
    return false;
}
