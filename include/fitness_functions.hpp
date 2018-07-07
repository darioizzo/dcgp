#ifndef DCGP_FITNESS_FUNCTIONS_H
#define DCGP_FITNESS_FUNCTIONS_H

#include <cmath>
#include <vector>

#include "expression.hpp"
#include <audi/functions.hpp>

namespace dcgp
{

/// Computes the quadratic error of a dCGP expression in approximating given data
template <typename T1, typename T2, typename T3>
T1 quadratic_error(const expression<T3> &ex, const std::vector<std::vector<T1>> &in_des,
                   const std::vector<std::vector<T2>> &out_des)
{
    using namespace std;

    T1 retval(0.);
    std::vector<T1> out_real;

    if (in_des.size() != out_des.size()) {
        throw std::invalid_argument("Size of the input vector must be the size of the output vector");
    }

    for (auto i = 0u; i < in_des.size(); ++i) {
        out_real = ex(in_des[i]);
        for (auto j = 0u; j < out_real.size(); ++j) {
            retval += (out_des[i][j] - out_real[j]) * (out_des[i][j] - out_real[j]);
        }
    }
    return retval / static_cast<int>(out_des.size());
}

} // namespace dcgp

#endif
