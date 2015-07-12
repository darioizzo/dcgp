#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional>

namespace dcgp {

using my_fun_type = std::function<double(double, double)>;
using my_d_fun_type = std::function<double(unsigned int, double, double)>;

struct basis_function
{
        basis_function(my_fun_type f, my_d_fun_type df) : m_f(f), m_df(df) {}

        double operator()(double x, double y) const
        {
                return m_f(x,y);
        }
        my_fun_type m_f;
        my_d_fun_type m_df;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
