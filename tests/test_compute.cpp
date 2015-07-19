#include <iostream>
#include <iomanip>
#include "../src/dcgp.h"

#define EPSILON 1e-13

bool f_fails(const dcgp::expression& p, const std::vector<double>& in, const std::vector<double>& out)
{
    std::vector<double> out_computed = p.compute(in);
    std ::cout << out_computed << std::endl;
    std ::cout << out << std::endl;
    for (auto i = 0u; i<out_computed.size(); ++i)
    {
        if (fabs(out_computed[i]-out[i]) > EPSILON) return true;
    }
    return false;
}

/// This test is passed when some predefined expressions encoded in predefined genes compute correctly. The data are hand-written.

int main() {
    /// Testing over Miller's test case from the PPSN 2014 tutorial
    dcgp::expression miller(2,4,2,3,4,dcgp::function_set::minimal);
    std::vector<unsigned int> x1({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});
    miller.set(x1);
    bool miller_test_fails = f_fails(miller, {2.,3.},{5,6,-18,0}) || f_fails(miller, {1.,-1.},{0,-1,-1,0}) || f_fails(miller, {-.123,2.345},{2.222,-0.288435,0.676380075,0});

    std::vector<unsigned int> x2({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 6, 5, 7, 3});
    miller.set(x2);
    miller_test_fails = miller_test_fails || f_fails(miller, {2.,3.},{-6,6,-18,0}) || f_fails(miller, {1.,-1.},{2,-1,-1,0}) || f_fails(miller, {-.123,2.345},{-4.69,-0.288435,0.676380075,0});

    /// Testing over a single row program
    dcgp::expression one_row(4,1,1,10,10,dcgp::function_set::minimal);
    x1 = {2, 3, 0, 0, 2, 2, 3, 0, 1, 1, 5, 4, 2, 6, 1, 0, 7, 7, 3, 6, 7, 1, 7, 6, 2, 4, 10, 2, 3, 2, 10};     ///(x/y)/(2z-(t*x))
    one_row.set(x1);
    bool one_raw_fails = f_fails(one_row, {2.,3.,4.,-2.}, {0.055555555555555552}) || f_fails(one_row, {-1.,1.,-1.,1.}, {1}) || f_fails(one_row,{0,1,2,3},{0});

    return miller_test_fails || one_raw_fails;
}



#undef EPSILON