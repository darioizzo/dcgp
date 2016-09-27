#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

#include <functional> // std::function
#include <utility>    // std::forward
#include <string>
#include <iostream>
#include <vector>
#include <audi/audi.hpp>

using namespace audi;
using gdual_d = audi::gdual<double>;

namespace dcgp {

/// Basis function
/**
 *
 *
 */
template<typename T>
struct basis_function
{
    /// Basic prototype of a kernel function for its simple evaluation
    using my_fun_type = std::function<T(const std::vector<T>&)>;
    /// Basic prototype of a kernel function for the evaluation of its printable form
    using my_print_fun_type = std::function<std::string(const std::vector<std::string>&)>;

    /// Constructor from std::function construction arguments
    /*
     * Construct a function that can be used as a kernel in a d-CGP expression
     *
     * @param[in] f constructs a dcgp::my_fun_type
     * @param[in] df constructs a dcgp::d_my_fun_type
     * @param[in] pf constructs a dcgp::my_print_fun_type
     * @param[in] name string containing the function name (ex. "sum")
     */
    template <typename U, typename V>
    basis_function(U &&f, V&&pf, std::string name):m_f(std::forward<U>(f)), m_pf(std::forward<V>(pf)), m_name(name) {}

    /// Parenthesis operator overload
    /**
    * Evaluates the basis_function in the point \p in
    *
    * @param[in] in evaluation point
    *
    * @return the function evaluated
    */
    T operator()(const std::vector<T>& in) const
    {
            return m_f(in);
    }
    T operator()(const std::initializer_list<T>& in) const
    {
            return m_f(in);
    }

    /// Parenthesis operator overload (std::string)
    /**
    * Returns a symbolic representation for the operation made by \f$f\f$
    *
    * @param[in] in vector of string with the symbolic names to be used
    *
    * @return the string representation of the operation (ex. "ln(s1+s2)")
    */
    std::string operator()(const std::vector<std::string>& in) const
    {
            return m_pf(in);
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
    /// Its symbolic representation
    my_print_fun_type m_pf;
    /// Its name
    std::string m_name;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
