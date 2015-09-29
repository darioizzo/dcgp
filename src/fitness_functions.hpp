#include <cmath>
#include <vector>

#include "expression.hpp"

namespace dcgp {

    /// Computes the error of a d-CGP expression in approximating given data
    template <typename T>
    T symbolic_regression(const expression& ex, 
        const std::vector<std::vector<T> >& in_des, 
        const std::vector<std::vector<T> >& out_des) 
    {
        using namespace std;

        T retval(0.);
        std::vector<T> out_real;

        if (in_des.size() != out_des.size())
        {
            throw std::invalid_argument("Size of the input vector must be the size of the output vector");
        }

        for (auto i = 0u; i < in_des.size(); ++i)
        {
            out_real = ex(in_des[i]);
            for (auto j = 0u; j < out_real.size(); ++j)
            {
                retval += abs(out_des[i][j] - out_real[j]);
            }
        }
        return retval;
    }
} // namespace
