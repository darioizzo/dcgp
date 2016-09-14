#include <iostream>

#include "../include/expression.hpp"
#include "../include/function_set.hpp"

// Here we search for first integrals of the Kepler's problem) using our "mutation suppression" method
// The hamiltonian is H = 1/(2m) (pr^2+pt^2 / r^2) + mu / r

double fitness(const dcgp::expression& ex, const std::vector<std::vector<double> >& in, double& check)
{
    double retval = 0;
    check = 0;
    for (auto i = 0u; i < in.size(); ++i) {
        auto T = ex.taylor(in[i], 1);               // We compute all the derivatives up to order one
        double dFpr= T[0].get_derivative({1, 0, 0, 0, 0, 0});
        double dFpt= T[0].get_derivative({0, 1, 0, 0, 0, 0});
        double dFqr= T[0].get_derivative({0, 0, 1, 0, 0, 0});
        double dFqt= T[0].get_derivative({0, 0, 0, 1, 0, 0});
        double pr = in[i][0];
        double pt = in[i][1];
        double qr = in[i][2];
        double qt = in[i][3];
        double m =  in[i][4];
        double mu = in[i][5];
        double err = dFpr * (pt*pt/m/qr/qr/qr - mu/qr/qr) + dFqr * (pr/m) + dFqt * (pt/qr/qr/m);
        retval += (err) * (err);                                 // We compute the quadratic error
        check += dFpr*dFpr + dFqr*dFqr + dFqt*dFqt; // If check is 0 then ex represent a constant or an ignorable variable
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
    dcgp::expression ex(6, 1, 1, 100, 50, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"pr","pt", "r", "th", "m", "mu"});

    // We create the grid over x
    std::vector<double> dumb(6);
    std::vector<std::vector<double> > in(50, dumb);
    for (auto i = 0u; i < in.size(); ++i) {
        in[i][0] = 0.12 + 20 / (in.size() - 1) * i;
        in[i][1] = 1 + 20 / (in.size() - 1) * i;
        in[i][2] = 0.10 + 20 / (in.size() - 1) * i;
        in[i][3] = 1 + 20 / (in.size() - 1) * i;
        in[i][4] = 0.11 + 0.9 / (in.size() - 1) * i;
        in[i][5] = 1 + 0.143 / (in.size() - 1) * i;
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
                    stream(std::cout, "New best found: gen: ", std::setw(7), gen, "\t value: ", newfits[i], "\n");
                    //std::cout << "Expression: " << ex(in_sym), "\n");
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
