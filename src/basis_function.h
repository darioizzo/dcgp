#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional>

namespace dcgp {

using my_fun_type = std::function<double(double, double)>;
using my_d_fun_type = std::function<double(unsigned int, double, double)>;
using my_print_fun_type = std::function<std::string(std::string, std::string)>;

struct basis_function
{

	// a slightly better version of basis_function(T f, U df):m_f(f),m_df(df) {} 
	// as it passes to the std::function constuctor exactly the same as received as inputs (i.e. no copies made)
	template <typename T, typename U, typename V>
    basis_function(T &&f, U &&df, V&&pf):m_f(std::forward<T>(f)),m_df(std::forward<U>(df)),m_pf(std::forward<V>(pf)) {}

    double operator()(double x, double y) const
    {
            return m_f(x,y);
    }
    std::string operator()(std::string x, std::string y) const
    {
            return m_pf(x,y);
    }
    my_fun_type m_f;
    my_d_fun_type m_df;
	my_print_fun_type m_pf;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
 