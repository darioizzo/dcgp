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
 * This class represents the function defining the generic CGP node. To be constructed
 * it accepts two functions having prototype \p T(const std::vector< \p T>&) and std::string(const std::vector<std::string>&)
 * computing, respectively the function value on generic inputs and the textual representation of the operation.
 *
 * The intended use would then be something like:
 * @code
 * basis_function<double> f(my_sum<double>, print_my_sum, "sum");
 * basis_function<double> f(my_sum<gdual_d>, print_my_sum, "sum");
 * @endcode
 *
 * tparam T The type of the function output (and inputs)
 */
template<typename T>
class basis_function
{
public:
    /// Basic prototype of a kernel function returning its evaluation
    using my_fun_type = std::function<T(const std::vector<T>&)>;
    /// Basic prototype of a kernel function returning its symbolic representation
    using my_print_fun_type = std::function<std::string(const std::vector<std::string>&)>;

    /// Constructor
    /*
     * Constructs a basis_function that can be used as kernel in a dCGP expression
     *
     * @param[in] f any callable with prototype T(const std::vector<T>&)
     * @param[in] pf any callable with prototype std::string(const std::vector<std::string>&)
     * @param[in] name string containing the function name (ex. "sum")
     *
     */
    template <typename U, typename V>
    basis_function(U &&f, V&&pf, std::string name):m_f(std::forward<U>(f)), m_pf(std::forward<V>(pf)), m_name(name) {}

    /// Parenthesis operator
    /**
    * Evaluates the basis_function in the point \p in
    *
    * @param[in] in the evaluation point as an std::vector<T>
    *
    * @return the function value
    */
    T operator()(const std::vector<T>& in) const
    {
            return m_f(in);
    }
    /// Parenthesis operator
    /**
    * Evaluates the basis_function in the point \p in
    *
    * @param[in] in the evaluation point as an std::initiaizer_list<T>
    *
    * @return the function value
    */
    T operator()(const std::initializer_list<T>& in) const
    {
            return m_f(in);
    }
    /// Parenthesis operator
    /**
    * Returns a symbolic representation of the operation made by \f$f\f$
    *
    * @param[in] in std::vector<std::string> with the symbolic names to be used
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
     * @param[in] d basis_function argument.
     *
     * @return reference to \p os.
     *
    */
    friend std::ostream& operator<<(std::ostream& os, const basis_function<T>& d)
    {
        os << d.m_name;
        return os;
    }

private:
    /// The function
    my_fun_type m_f;
    /// Its symbolic representation
    my_print_fun_type m_pf;
    /// Its name
    std::string m_name;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
