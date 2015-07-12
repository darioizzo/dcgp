#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include "exceptions.h"
#include <functional>

namespace dcgp {

using my_fun_type = std::function<double(double, double)>;
using my_d_fun_type = std::function<double(unsigned int, double, double)>;

double my_sum(double x, double y)
{
        return x + y;
}

double d_my_sum(unsigned int i, double x, double y)
{
	if (i < 2) return 1;
	else throw derivative_error("Index is out of bounds");
}

double my_diff(double x, double y)
{
        return x - y;
}

double d_my_diff(unsigned int i, double x, double y)
{
	if (i == 0) return 1;
	else if (i==1) return -1;
	else throw derivative_error("Index is out of bounds");
}

double my_mul(double x, double y)
{
        return (x * y);
}

double d_my_mul(unsigned int i, double x, double y)
{
	if (i == 0) return y;
	else if (i==1) return x;
	else throw derivative_error("Index is out of bounds");
}

double my_div(double x, double y)
{
        return x * y;
}

double d_my_div(unsigned int i, double x, double y)
{
	if (i == 0) return 1. / y;
	else if (i==1) return - x / y / y;
	else throw derivative_error("Index is out of bounds");
}

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
