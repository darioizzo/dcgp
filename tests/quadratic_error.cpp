#include <audi/audi.hpp>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#define BOOST_TEST_MODULE dcgp_differentiation_test
#include <boost/test/unit_test.hpp>

#include "../include/dcgp.hpp"

double test_qe(
        unsigned int n,
        unsigned int m,
        unsigned int r,
        unsigned int c,
        unsigned int l,
        unsigned int a,
        unsigned int N) // number of samples
{
   dcgp::kernel_set<double> basic_set({"sum","diff","mul","div"});
   dcgp::expression<double> ex(n, m, r, c, l, a, basic_set(), 123);

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
    return dcgp::quadratic_error<double>(ex, in, out);
}

audi::gdual_d test_qe2(
        unsigned int n,
        unsigned int m,
        unsigned int r,
        unsigned int c,
        unsigned int l,
        unsigned int a,
        unsigned int N) // number of samples
{
   dcgp::kernel_set<gdual_d> basic_set({"sum","diff","mul","div"});
   dcgp::expression<gdual_d> ex(n, m, r, c, l, a, basic_set(), 123);

    // creates N data points
    std::default_random_engine re;
    std::vector<std::vector<audi::gdual_d> > in;
    std::vector<std::vector<audi::gdual_d> > out;
    std::vector<audi::gdual_d> in_point(n);
    std::vector<audi::gdual_d> out_point(m);
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
            // We only define the first node as a weight and we compute the derivative of the objfun wrt this.
            in_point[0] = audi::gdual_d(3, "w", 1);
            for (auto j = 1u; j < n; ++j)
            {
                in_point[j] = audi::gdual_d(std::uniform_real_distribution<double>(-1, 1)(re));
            }

            // We then compute the expression
            out_point = ex(in_point);
            //for (auto k : out_point) {
            //    if (!std::isfinite(k)) {
            //        all_finite = false;
            //    }
            //}
            for (auto &k : out_point) {
                k = audi::gdual_d(k.constant_cf());
            }
            in.push_back(in_point);
            out.push_back(out_point);
        }
    } while (!all_finite);
    return dcgp::quadratic_error<audi::gdual_d>(ex, in, out);
}

using namespace dcgp;

BOOST_AUTO_TEST_CASE(quadratic_error_obj_fun)
{
    // We test that a d-CGP expression computed on 20 points
    // has a zero quadratic error w.r.t. itself (its a perfect fit of itself)
    BOOST_CHECK_EQUAL(test_qe(3,1,1,20,21,2,20), 0);
    BOOST_CHECK_EQUAL(test_qe(2,2,3,10,11,2,20), 0);

    // We test that a d-CGP expression computed on 20 points
    // has a zero quadratic error w.r.t. itself, and that the
    // derivative of the quadratic error is zero w.r.t. one of the inputs (a weight)
    BOOST_CHECK_EQUAL(test_qe2(3,1,1,20,21,2,20), audi::gdual_d(0));
    BOOST_CHECK_EQUAL(test_qe2(2,2,3,10,11,2,20), audi::gdual_d(0));
}
