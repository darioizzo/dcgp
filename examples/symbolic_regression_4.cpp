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

// In this second tutorial we show how to find an exact formula for some input data where we also
// want to learn some constants.
//
// This case is more generic and interesting than the one treated in the previous tutorial,
// and can be considered, already, of use in industrial settings.
//
// NOTE: constants can, in general, be learned via two main techniques
// 1) evolutionary (common practice in GP)
// 2) 'backpropagation' (this is original with dCGP)
//
// In this tutorial we use 2)
//
// We use the problem P2 from the dcgp::gym, that is x**5 - pi*x**3 + 2 pi / x

int main()
{
    // We load the data (using problem P1 from the gym)
    std::vector<std::vector<double>> X, Y;
    gym::generate_vladi6(X, Y);

    // We instantiate a symbolic regression problem with one only ephemeral constants.
    auto n_eph = 1u;
    symbolic_regression udp(X, Y, 1u, 100, 101, 2u, kernel_set<double>({"sum", "diff", "mul", "sin", "cos"})(), n_eph, true);

    // We init a population with four individuals.
    pagmo::population pop{udp, 100};

    // We instantiate the memetic solver setting 1000 maximum generation and one active mutation (minimum)
    dcgp::momes4cgp uda{100, 4};
    pagmo::algorithm algo{uda};
    algo.set_verbosity(1u);

    // We solve
    pop = algo.evolve(pop);
    pagmo::print("\n");
    // We print on screen the non dominated front
    auto ndf = pagmo::non_dominated_front_2d(pop.get_f());
    for (decltype(ndf.size()) i = 0u; i < ndf.size(); ++i) {
        auto idx = ndf[i];
        auto prettier = udp.prettier(pop.get_x()[idx]);
        trim_left_if(prettier, is_any_of("["));
        trim_right_if(prettier, is_any_of("]"));
        pagmo::print("Loss: ", pop.get_f()[idx][0], std::setw(15), "Complexity: ", pop.get_f()[idx][1], std::setw(10), "Formula:", prettier,"\n");
    }

    return false;
}
