#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional>

namespace dcgp {

using my_fun_type = std::function<double(double, double)>;
using my_d_fun_type = std::function<double(unsigned int, double, double)>;

struct basis_function
{

	// a slightly better version of basis_function(T f, U df):m_f(f),m_df(df) {} 
	// as it passes to the std::function constuctor exactly the same as received as inputs (i.e. no copies made)
	template <typename T, typename U>
    basis_function(T &&f, U &&df):m_f(std::forward<T>(f)),m_df(std::forward<U>(df)) {}

    double operator()(double x, double y) const
    {
            return m_f(x,y);
    }
    my_fun_type m_f;
    my_d_fun_type m_df;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
 