#include <iostream>

#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>

// Here we search for first integrals of the Kepler's problem) using our "mutation suppression" method
// The hamiltonian is H = 1/(2m) (pr^2+pt^2 / r^2) + mu / r

double fitness(const dcgp::expression<gdual_d> &ex, const std::vector<std::vector<gdual_d>> &in, double &check)
{
    double retval = 0;
    check = 0;
    for (auto i = 0u; i < in.size(); ++i) {
        auto T = ex(in[i]); // We compute all the derivatives up to order one
        std::vector<std::string> symbol_list{"pr", "pt", "r", "th", "m", "mu"};
        for (auto sym : symbol_list) {
            T[0] += gdual_d(0, sym, 0); // We make sure that the symbols are all in the final expression
        }
        double dFpr = T[0].get_derivative({{"dpr", 1}});
        // double dFpt= T[0].get_derivative({{"dpt", 1}});
        double dFqr = T[0].get_derivative({{"dr", 1}});
        double dFqt = T[0].get_derivative({{"dth", 1}});
        double pr = in[i][0].constant_cf();
        double pt = in[i][1].constant_cf();
        double qr = in[i][2].constant_cf();
        // double qt = in[i][3].constant_cf();
        double m = in[i][4].constant_cf();
        double mu = in[i][5].constant_cf();
        double err = dFpr * (pt * pt / m / qr / qr / qr - mu / qr / qr) + dFqr * (pr / m) + dFqt * (pt / qr / qr / m);
        retval += (err) * (err); // We compute the quadratic error
        check += dFpr * dFpr + dFqr * dFqr
                 + dFqt * dFqt; // If check is 0 then ex represent a constant or an ignorable variable
    }

    return retval;
}

using namespace dcgp;

int main()
{
    // Random seed
    std::random_device rd;

    // Function set
    dcgp::kernel_set<gdual_d> basic_set({"sum", "diff", "mul", "div"});

    // d-CGP expression
    dcgp::expression<gdual_d> ex(6, 1, 1, 100, 50, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"pr", "pt", "r", "th", "m", "mu"});

    // We create the grid over x
    std::vector<std::vector<gdual_d>> in(50u);
    for (auto i = 0u; i < in.size(); ++i) {
        gdual_d pr_var(0.12 + 20. / static_cast<double>((in.size() - 1)) * i, "pr", 1u);
        gdual_d pt_var(1. + 20. / static_cast<double>((in.size() - 1)) * i, "pt", 1u);
        gdual_d r_var(0.12 + 20. / static_cast<double>((in.size() - 1)) * i, "r", 1u);
        gdual_d th_var(1. + 20. / static_cast<double>((in.size() - 1)) * i, "th", 1u);
        gdual_d m_var(0.12 + 20. / static_cast<double>((in.size() - 1)) * i, "m", 1u);
        gdual_d mu_var(1. + 20. / static_cast<double>((in.size() - 1)) * i, "mu", 1u);
        in[i] = std::vector<gdual_d>{pr_var, pt_var, r_var, th_var, m_var, mu_var};
    }

    // We run the (1-4)-ES
    double best_fit = 1e32;
    std::vector<double> newfits(4, 0.);
    std::vector<std::vector<unsigned int>> newchromosomes(4);
    std::vector<unsigned int> best_chromosome(ex.get());
    unsigned int gen = 0;
    do {
        gen++;
        for (auto i = 0u; i < newfits.size(); ++i) {
            ex.set(best_chromosome);
            double check = 0;
            while (check < 1e-3) {
                ex.mutate_active(6);
                newfits[i] = fitness(ex, in, check); // Total fitness
            }
            newchromosomes[i] = ex.get();
        }
        for (auto i = 0u; i < newfits.size(); ++i) {
            if (newfits[i] <= best_fit) {
                if (newfits[i] != best_fit) {
                    stream(std::cout, "New best found: gen: ", std::setw(7), gen, "\t value: ", newfits[i], "\n");
                    // std::cout << "Expression: " << ex(in_sym), "\n");
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
    stream(std::cout, "Point: ", in[2], "\n");
    stream(std::cout, "Taylor: ", ex(in[2]), "\n");
}
