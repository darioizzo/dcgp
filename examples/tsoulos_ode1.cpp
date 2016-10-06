#include <iostream>

#include "../include/expression.hpp"
#include "../include/function_set.hpp"

// Here we solve the differential equation dy = (2x - y) / x from Tsoulos paper
// Tsoulos and Lagaris: "Solving Differential equations with genetic programming"

double fitness(const dcgp::expression<gdual_d>& ex, const std::vector<std::vector<gdual_d> >& in)
{
    double retval = 0;
    for (auto i = 0u; i < in.size(); ++i) {
        auto T = ex(in[i]);                   // We compute the expression and thus the derivatives
        double y = T[0].get_derivative({0});
        double dy = T[0].get_derivative({1});
        double x = in[i][0].constant_cf();
        double ode1 = (2. * x - y) / x;
        retval += (ode1 - dy) * (ode1 - dy);            // We compute the quadratic error
    }
    return retval;
}

int main () {
    // Random seed
    std::random_device rd;

    // Function set
    dcgp::function_set<gdual_d> basic_set({"sum", "diff", "mul", "div"});

    // d-CGP expression
    dcgp::expression<gdual_d> ex(1, 1, 1, 15, 16, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"x"});

    // We create the grid over x
    std::vector<std::vector<gdual_d> > in(10u);
    for (auto i = 0u; i < in.size(); ++i) {
        gdual_d point(0.1 + 0.9 / (in.size() - 1) * i, "x", 1);
        in[i].push_back(point);   // 1, .., 2
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
            ex.mutate_active(2);
            auto fitness_ic = ex(std::vector<gdual_d>{gdual_d(1.)})[0] - 3.;  // Penalty term to enforce the initial conditions
            newfits[i] = fitness(ex, in) + fitness_ic.constant_cf() * fitness_ic.constant_cf();     // Total fitness
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
    } while (best_fit > 1e-3 && gen < 10000);
    dcgp::stream(std::cout, "Number of generations: ", gen, "\n");
    dcgp::stream(std::cout, "Expression: ", ex(in_sym), "\n");
}
