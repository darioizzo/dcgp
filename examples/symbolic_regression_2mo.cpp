#include <boost/algorithm/string.hpp>
#include <chrono>
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

#include <dcgp/algorithms/moes4cgp.hpp>
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
    symbolic_regression udp(X, Y, 1u, 20u, 21u, 2u, kernel_set<double>({"sum", "diff", "mul", "pdiv"})(), n_eph, true,
                            0u, "MSE");

    // We init a population with four individuals
    pagmo::population pop{udp, 4u};

    // We use an ES startegy also to learn constants
    dcgp::moes4cgp uda(100000u, 4u, 1e-8, true);
    pagmo::algorithm algo{uda};
    algo.set_verbosity(500u);

    auto t1 = std::chrono::high_resolution_clock::now();
    pop = algo.evolve(pop);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    std::cout << duration;

    // Finally we print on screen the non dominated front.
    pagmo::print("\nNon dominated Front at the end:\n");
    auto ndf = pagmo::non_dominated_front_2d(pop.get_f());
    for (decltype(ndf.size()) i = 0u; i < ndf.size(); ++i) {
        auto idx = ndf[i];
        auto prettier = udp.prettier(pop.get_x()[idx]);
        trim_left_if(prettier, is_any_of("["));
        trim_right_if(prettier, is_any_of("]"));
        pagmo::print(std::setw(2), i + 1, " - Loss: ", std::setw(13), std::left, pop.get_f()[idx][0], std::setw(15),
                     "Complexity: ", std::left, std::setw(5), pop.get_f()[idx][1], std::setw(10), "Formula: ", prettier,
                     "\n");
    }
    return false;
}
