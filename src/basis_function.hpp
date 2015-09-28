#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional> // std::function
#include <utility>    // std::forward
#include <string>
#include <iostream>
#include <vector>
#include <audi/gdual.hpp>

namespace dcgp {

/// Basic prototype of a kernel function for its simple evaluation
using my_fun_type = std::function<double(const double&, const double&)>;
/// Basic prototype of a kernel function for the evaluation of its Taylor expansion
using d_my_fun_type = std::function<gdual(const gdual&, const gdual&)>;
/// Basic prototype of a kernel function for the evaluation of its printable form
using my_print_fun_type = std::function<std::string(std::string, std::string)>;

/// Basis function
/**
 * This struct represent one of the kernel functions to be used in a d-CGP expresion. 
 * It contains three std::function, whose type are my_fun_type, my_d_fun_type,
 * my_print_fun_type. These allow to compute the function value, its Taylor expansion
 * and its symbolic representation.
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 * @author Francesco Biscani (bluescarni@gmail.com)
 */
struct basis_function
{
    /// Constructor from std::function construction arguments
    /*   
     * Construct a function that can be used as a kernel in a d-CGP expression
     *
     * @param[in] f constructs a dcgp::my_fun_type 
     * @param[in] df constructs a dcgp::d_my_fun_type 
     * @param[in] pf constructs a dcgp::my_print_fun_type 
     * @param[in] name string containing the function name (ex. "sum")
     */
    template <typename T, typename U, typename V>
    basis_function(T &&f, U &&df, V&&pf, std::string name):m_f(std::forward<T>(f)), m_df(std::forward<U>(df)), m_pf(std::forward<V>(pf)), m_name(name) {}

    /// Parenthesis operator overload (double)
    /**
    * Evaluates \f$f(x, y)\f$, that is the basis_function in a point \p x, \p y
    *
    * @param[in] x first input
    * @param[in] y second input
    *
    * @return the function evaluated in \f$x,y\f$
    */
    double operator()(double x, double y) const
    {
            return m_f(x,y);
    }

    /// Parenthesis operator overload (audi::gdual)
    /**
    * Computes, in the algebra of truncated polynomial, the function \f$f\f$. The
    * result will be the Taylor expansion of the function \f$f\f$ composed with
    * the functions \p p1 and \p p2
    *
    * @param[in] p1 first input
    * @param[in] p2 second input
    *
    * @return the Taylor representation of \f$f\f$
    */
    audi::gdual operator()(const audi::gdual & p1, const audi::gdual & p2) const
    {
            return m_df(p1,p2);
    }

    /// Parenthesis operator overload (std::string)
    /**
    * Returns a symbolic representation for the operation made by \f$f\f$
    *
    * @param[in] s1 first input
    * @param[in] s2 second input
    *
    * @return the string representation of the operation (ex. "ln(s1+s2)")
    */
    std::string operator()(std::string s1, std::string s2) const
    {
            return m_pf(s1,s2);
    }

    /// Overloaded stream operator
    /**
     * Will stream the function name
     * 
     * @param[in,out] os target stream.
     * @param[in] d dcgp::basis_function argument.
     * 
     * @return reference to \p os.
     * 
    */
    friend std::ostream& operator<<(std::ostream& os, const basis_function& d)
    {
        os << d.m_name;
        return os;
    }

    /// The function
    my_fun_type m_f;
    /// Its derivatives
    d_my_fun_type m_df;
    /// Its symbolic representation
    my_print_fun_type m_pf;
    /// Its name
    std::string m_name;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
 