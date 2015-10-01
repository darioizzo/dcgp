#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#include "../src/expression.hpp"
#include "../src/function_set.hpp"
#include "../src/expression.hpp"
#include "../src/fitness_functions.hpp"

bool kq(
        unsigned int r,
        unsigned int c,
        unsigned int l,
        unsigned int N) // number of samples
{
    // Random seed
    std::random_device rd;

    // Function set
    dcgp::function_set basic_set({"sum", "diff", "mul", "div"});

    // d-CGP expression
    dcgp::expression ex(1, 1, r, c, l, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"x"});

    // 1) we create N data points 
    std::default_random_engine re(12);
    std::vector<std::vector<double> > in;
    std::vector<std::vector<double> > out;
    std::vector<double> in_point(1);
    std::vector<double> out_point(1);

    // Koza quintic polynomial x^5 - 2x^3 + x
    std::function<double(double)> f = [](double x){return x*x*x*x*x - 2*x*x*x + x;};
    for (auto i = 0u; i < N; ++i)
    {
        in_point[0] = std::uniform_real_distribution<double>(-1, 1)(re);
        out_point[0] = f(in_point[0]);
        in.push_back(in_point);
        out.push_back(out_point);
    }

    /// 2) we use a simple ES(1+4) to evolve an expression that represents our target
    /// note that the problem is a maximization problem
    double best_fit = 1e32;
    std::vector<double> newfits(4,0.);
    std::vector<std::vector<unsigned int> > newchromosomes(4);
    std::vector<unsigned int> best_chromosome;
    unsigned int gen = 0;

    do
    {
        gen++;
        for (auto i = 0u; i < newfits.size(); ++i) {
            ex.mutate_active();
            newfits[i] = quadratic_error(ex, in, out);
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
    } while (best_fit > 1e-3);
    ex.set(best_chromosome);

    std::cout << "Number of generations: " << gen << std::endl;
    std::cout << "Final expression: " << ex(in_sym) << std::endl;
    return false;
}


int main() {
    return kq(1,50,51,10);
}

#undef TOL
