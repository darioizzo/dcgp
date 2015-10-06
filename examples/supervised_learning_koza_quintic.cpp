#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#include "../src/function_set.hpp"
#include "../src/expression.hpp"
#include "../src/fitness_functions.hpp"
#include "detail/es.hpp"

bool kq(
        unsigned int r,
        unsigned int c,
        unsigned int l,
        unsigned int a,
        unsigned int N) // number of samples
{
    // Random seed
    std::random_device rd;

    // Function set
    dcgp::function_set basic_set({"sum", "diff", "mul", "div"});

    // d-CGP expression
    dcgp::expression ex(1, 1, r, c, l, a, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"x"});

    // 1) we create N data points for the Koza quintic polynomial x^5 - 2x^3 + x
    std::default_random_engine re(12);
    std::vector<std::vector<double> > in;
    std::vector<std::vector<double> > out;
    std::vector<double> in_point(1);
    std::vector<double> out_point(1);
    std::function<double(double)> f = [](double x){return x*x*x*x*x - 2*x*x*x + x;};
    for (auto i = 0u; i < N; ++i)
    {
        in_point[0] = std::uniform_real_distribution<double>(-1, 1)(re);
        out_point[0] = f(in_point[0]);
        in.push_back(in_point);
        out.push_back(out_point);
    }

    // 2) we use a simple ES(1+4) to evolve an expression that represents our target. 
    // Mutation only mutates 2 active genes
    es_params params{4, "active", 2, 0};
    es(in, out, ex, params);

    std::cout << "Final expression: " << ex(in_sym) << std::endl;
    return false;
}


int main() {
    return kq(1,15,16,2,10);
}

