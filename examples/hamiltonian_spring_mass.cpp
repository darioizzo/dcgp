#include <iostream>

#include "../src/expression.hpp"
#include "../src/function_set.hpp"

// Here we search for first integrals of the mass spring sistem (one dimension) using our "mutation suppression" method
// The hamiltonian is H = 1/2 p^2 + 1/2 q^2

double fitness(const dcgp::expression& ex, const std::vector<std::vector<double> >& in, double& check)
{
    double retval = 0;
    check = 0;
    for (auto i = 0u; i < in.size(); ++i) {
        auto T = ex.taylor(in[i], 1);                   // We compute all the derivatives up to order one
        double dFp = T[0].get_derivative({1, 0});
        double dFq = T[0].get_derivative({0, 1});
        double p = in[i][0];
        double q = in[i][1];
        double err = - dFp * q + dFq * p;                       
        retval += (err) * (err);  
        check += dFp*dFp + dFq*dFq;                      // We compute the quadratic error 
    }

    return retval;
}

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
            double check = 0;
            while (check < 1e-3) {
                ex.mutate_active(6);
                newfits[i] = fitness(ex, in, check);   // Total fitness
            }
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

    std::cout << "Number of generations: " << gen << std::endl;
    std::cout << "Expression: " <<  ex << std::endl;
    std::cout << "Expression: " << ex(in_sym) << std::endl;
    std::cout << "Point: " <<  in[2] << std::endl;
    std::cout << "Taylor: " <<  ex.taylor(in[2], 1) << std::endl;

}