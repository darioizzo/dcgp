#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#define BOOST_TEST_MODULE dcgp_differentiation_test
#include <boost/test/unit_test.hpp>

#include "../src/dcgp.hpp"

double test_sr(
        unsigned int n,
        unsigned int m,
        unsigned int r,
        unsigned int c,
        unsigned int l,
        unsigned int N) // number of samples
{
   dcgp::function_set basic_set({"sum","diff","mul","div"});
   dcgp::expression ex(n, m, r, c, l, basic_set(), 123);

    // creates N data points 
    std::default_random_engine re;
    std::vector<std::vector<double> > in;
    std::vector<std::vector<double> > out;
    std::vector<double> in_point(n);
    std::vector<double> out_point(m);
    bool all_finite;
    // assuming that mutating the program will sooner or later produce an expression that is not nan or inf in any of the points ...
    // ... the following loop is not infinite
    do 
    {
        ex.mutate_active();
        all_finite=true;
        in.clear();
        out.clear();
        for (auto i = 0u; i < N; ++i)
        {
            for (auto j = 0u; j < n; ++j) 
            {
                in_point[j] = std::uniform_real_distribution<double>(-1, 1)(re);
            }
            out_point = ex(in_point);
            for (auto k : out_point) {
                if (!std::isfinite(k)) {
                    all_finite = false;
                }
            }
            in.push_back(in_point);
            out.push_back(out_point);
        }
    } while (!all_finite);
    return dcgp::symbolic_regression<double>(ex, in, out);
}

using namespace dcgp;

BOOST_AUTO_TEST_CASE(symbolic_regression_obj_fun)
{
    double res = test_sr(3,1,1,20,21,100);
    std::cout << res << std::endl;
}