#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <string>
#include "../src/dcgp.h"

std::vector<double> numeric_d(const dcgp::expression& ex, unsigned int wrt,  const std::vector<double>& in)
{
    if (ex.get_n() != in.size()) throw dcgp::input_error("The dimensions of the input point vector is not equal to the d-CGP inputs number");
    if (wrt >= ex.get_n()) throw dcgp::input_error("Trying to make a derivative with respect to a non existing variable");
    std::vector<double> retval(ex.get_m());
    double h = 1e-8;
    std::vector<double> in_plus_h = in;
    std::vector<double> in_minus_h = in;
    in_plus_h[wrt] += h;
    in_minus_h[wrt] -= h;
    std::vector<double> out_plus_h = ex(in_plus_h);
    std::vector<double> out_minus_h = ex(in_minus_h);
    for (auto i = 0u; i < retval.size(); ++i)
    {
        retval[i] = (out_plus_h[i] - out_minus_h[i]) / 2. / h;
    }
    return retval;
}

std::vector<double> numeric_d2(const dcgp::expression& ex, unsigned int wrt,  const std::vector<double>& in)
{
    if (ex.get_n() != in.size()) throw dcgp::input_error("The dimensions of the input point vector is not equal to the d-CGP inputs number");
    if (wrt >= ex.get_n()) throw dcgp::input_error("Trying to make a derivative with respect to a non existing variable");
    std::vector<double> retval(ex.get_m());
    double h = 1e-4; // numerical differentiation does not allow for smaller steps -> error becomes huge
    std::vector<double> in_plus_h = in;
    std::vector<double> in_minus_h = in;
    in_plus_h[wrt] += h;
    in_minus_h[wrt] -= h;
    std::vector<double> out = ex(in);
    std::vector<double> out_plus_h = ex(in_plus_h);
    std::vector<double> out_minus_h = ex(in_minus_h);

    for (auto i = 0u; i < retval.size(); ++i) // http://www2.math.umd.edu/~dlevy/classes/amsc466/lecture-notes/differentiation-chap.pdf
    {
        retval[i] = (out_plus_h[i] - 2 * out[i] + out_minus_h[i]) / h / h;
    }
    return retval;
}

bool test_fails(
        unsigned int n,
        unsigned int m,
        unsigned int r,
        unsigned int c,
        unsigned int l,
        std::vector<dcgp::basis_function> f_set,
        unsigned int N
        ) // number of samples
{
    // A random expression
    dcgp::expression ex(n, m, r, c, l, f_set);
    // We create the symbolic variables in case output is needed
    std::vector<std::string> in_sym;
    for (auto i = 0u; i<n; ++i)
    {
        in_sym.push_back("x" + std::to_string(i));
    }
    for (auto trial = 0u; trial < N;++trial) {
        // Pick a random input value in [-1, 1)
        std::vector<double> random_in(ex.get_n());
        std::default_random_engine re;
        for (auto i = 0u; i < random_in.size(); ++i)
        {
            random_in[i] = std::uniform_real_distribution<double>(-1, 1)(re);
        }

        // For each of the input variables
        for (auto i = 0u; i < n; ++i)
        {
            // Compute f, f' and f'' using automated differentiation
            std::vector<std::vector<double> > jet = ex.compute_d2(i,2,random_in);
            // Compute f' anf f'' using numerical differentiation
            std::vector<double> df = numeric_d(ex, i, random_in);
            std::vector<double> ddf = numeric_d2 (ex, i, random_in);
            // Compare the results
            for (auto k = 0u; k < df.size(); ++k)
            {
                if (fabs(jet[1][k] - df[k]) > std::max(1e-3 * fabs(jet[1][k]), 1e-3)) {
                    std::cout << "Failed First derivative: " << "\n";
                    std::cout << "Expression: " << "\n";
                    std::cout << ex(in_sym) << "\n";
                    std::cout << "Automated differentiation: " << jet[1] << "\n";
                    std::cout << "Numerical differentiation: " << df << "\n";
                    return true;
                }
                if (fabs(jet[2][k] - ddf[k]) > std::max(1e-3 * fabs(jet[2][k]), 1e-3)) {    
                    std::cout << "Failed Second derivative: " << "\n";
                    std::cout << "Expression: " << "\n";
                    std::cout << ex(in_sym) << "\n";
                    std::cout << "Automated differentiation: " << jet[2] << "\n";
                    std::cout << "Numerical differentiation: " << ddf << "\n";
                    return true;
                }
            }
        }
    }
   return false;
}

/// This test compares numerical and automated differentiation and passes if they compare well to tolerance
/// Note that the tolerance is set to be rather big as numerical differentiation sucks big time
int main() {
    std::vector<dcgp::basis_function> function_set = dcgp::function_set::minimal;
    return test_fails(2,4,2,3,4, function_set, 1000) ||
           test_fails(2,4,10,10,11, function_set, 1000) ||
           test_fails(2,4,20,20,21, function_set, 1000) ||
           test_fails(1,1,1,100,101, function_set, 1000) ||
           test_fails(1,1,2,100,101, function_set, 1000) ||
           test_fails(1,1,3,100,101, function_set, 1000);
}

