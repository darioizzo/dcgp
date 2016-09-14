#include <iostream>

#include "../include/expression.hpp"
#include "../include/function_set.hpp"

// Here we search for first integrals of a mass spring sistem (one dimension) using Lipson method
// The hamiltonian is H = p^2 + q^2 and is consistently found by the evolution

double fitness(const dcgp::expression& ex, const std::vector<std::vector<double> >& in)
{
    double retval = 0;
    for (auto i = 0u; i < in.size(); ++i) {
        auto T = ex.taylor(in[i], 1);                   // We compute all the derivatives up to order one
        double dFp = T[0].get_derivative({1, 0});
        double dFq = T[0].get_derivative({0, 1});
        double p = in[i][0];
        double q = in[i][1];
        double err = dFp / dFq - p / q;                 // Here we set (dp/dt) / (dq/dt) = dp/dq
        retval += (err) * (err);
    }

    return retval;
}

using namespace dcgp;

int main () {
    // Random seed
    std::random_device rd;

    // Function set
    dcgp::function_set basic_set({"sum", "diff", "mul", "div"});

    // d-CGP expression
    dcgp::expression ex(2, 1, 1, 15, 16, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"p","q"});

    // We create the grid over x
    std::vector<double> dumb(2);
    std::vector<std::vector<double> > in(10, dumb);
    for (auto i = 0u; i < in.size(); ++i) {
        in[i][0] = 0.12 + 0.9 / (in.size() - 1) * i; // 0.1, .., 1
        in[i][1] = 1 - 0.143 / (in.size() - 1) * i; // 1, 0.9, .. , 0.1
    }
    // We run the (1-4)-ES
    double best_fit = 1e32;
    std::vector<double> newfits(4, 0.);
    std::vector<std::vector<unsigned int> > newchromosomes(4);
    std::vector<unsigned int> best_chromosome(ex.get());
    unsigned int gen = 0;
    do
    {
        gen++;
        for (auto i = 0u; i < newfits.size(); ++i) {
            ex.set(best_chromosome);
            ex.mutate_active(6);
            newfits[i] = fitness(ex, in);   // Total fitness
            newchromosomes[i] = ex.get();
        }
        for (auto i = 0u; i < newfits.size(); ++i) {
            if (newfits[i] <= best_fit) {
                if (newfits[i] != best_fit) {
                    std::cout << "New best found: gen: " << std::setw(7) << gen << "\t value: " << newfits[i] << std::endl;
                    //std::cout << "Expression: " << ex(in_sym) << std::endl;
                }
                best_fit = newfits[i];
                best_chromosome = newchromosomes[i];
                ex.set(best_chromosome);
            }
        }
    } while (best_fit > 1e-12 && gen < 10000);

    stream(std::cout, "Number of generations: ", gen, "\n");
    stream(std::cout, "Expression: ", ex, "\n");
    stream(std::cout, "Expression: ", ex(in_sym), "\n");
    stream(std::cout, "Point: ",  in[2], "\n");
    stream(std::cout, "Taylor: ",  ex.taylor(in[2], 1), "\n");

}
