#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional>
#include <string>

namespace dcgp {

using my_fun_type = std::function<double(double, double)>;
using my_d_fun_type = std::function<double(unsigned int, double, double)>;
using my_print_fun_type = std::function<std::string(std::string, std::string)>;

/// Basis function
/**
 * This struct represent a generic function in d-CGP. It contains the std::function, whose type is
 * defined in my_fun_type, my_d_fun_type, my_print_fun_type, that allow to compute the function
 * value, its derivatives and its symbolic representation
 *
 * All functions that are represented in a d-CGP encoding must derive from this class and thus
 * the implementation of the function, its derivative and its symbolic representation must be
 * available as dcgp::my_fun_type, dcgp::my_d_fun_type and dcgp::my_print_fun_type in order
 * to be able to construct this object
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
struct basis_function
{
    /// Constructor from std::function construction arguments
    template <typename T, typename U, typename V>
    basis_function(T &&f, U &&df, V&&pf):m_f(std::forward<T>(f)),m_df(std::forward<U>(df)),m_pf(std::forward<V>(pf)) {}

    /// Overload of operator(double, double)
    /**
    * Allows to call a dcgp::basis_function with the syntax f(double x, double y) and get
    * the function value in x,y in return
    */
    double operator()(double x, double y) const
    {
            return m_f(x,y);
    }

    /// Overload of operator(std::string, std::string)
    /**
    * Allows to call a dcgp::basis_function with the syntax f(std::string x, std::string y) and get
    * a symbolic representation of the function in the variables std::string x, std::string y
    */
    std::string operator()(std::string x, std::string y) const
    {
            return m_pf(x,y);
    }

    /// The function
    my_fun_type m_f;
    /// Its derivatives
    my_d_fun_type m_df;
    /// Its symbolic representation
    my_print_fun_type m_pf;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
 