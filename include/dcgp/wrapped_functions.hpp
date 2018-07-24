#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <audi/audi.hpp>
#include <audi/functions.hpp>
#include <string>
#include <vector>

#include <dcgp/type_traits.hpp>

using namespace audi;

namespace dcgp
{

// SFINAE dust (to hide under the carpet). Its used to enable the templated
// version of the various functions that can construct a kernel object. Only for
// double and a gdual type Complex could also be allowed.
template <typename T>
using f_enabler = typename std::enable_if<std::is_same<T, double>::value || is_gdual<T>::value, int>::type;

// Allows to overload in templates std functions with audi functions
using namespace audi;

/*--------------------------------------------------------------------------
 *                                  N-ARITY FUNCTIONS
 *------------------------------------------------------------------------**/
template <typename T, f_enabler<T> = 0>
T my_sum(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return retval;
}

std::string print_my_sum(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
T my_diff(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval -= in[i];
    }
    return retval;
}

std::string print_my_diff(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "-" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
T my_mul(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval *= in[i];
    }
    return retval;
}

std::string print_my_mul(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "*" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
T my_div(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval /= in[i];
    }
    return retval;
}

std::string print_my_div(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "/" + in[i];
    }
    return "(" + retval + ")";
}

/*--------------------------------------------------------------------------
 *                                  Suitable for ANN
 *------------------------------------------------------------------------**/

// sigmoid function: 1 / (1 + exp(- (a + b + c + d+ .. + ))
template <typename T, f_enabler<T> = 0>
T my_sig(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return 1. / (1. + audi::exp(-retval));
}

std::string print_my_sig(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "sig(" + retval + ")";
}

// tanh function:
template <typename T, f_enabler<T> = 0>
T my_tanh(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return audi::tanh(retval);
}

std::string print_my_tanh(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "tanh(" + retval + ")";
}

// ReLu function:
template <typename T, f_enabler<T> = 0>
T my_relu(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    (retval < 0) ? retval= 0 : retval=retval ;
    return retval;
}

std::string print_my_relu(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "ReLu(" + retval + ")";
}

/*--------------------------------------------------------------------------
 *                                  BINARY FUNCTIONS
 *------------------------------------------------------------------------**/

// protected division: returns 1 if the the dividend equals the divisor
template <typename T, f_enabler<T> = 0>
T my_pdiv(const std::vector<T> &in)
{
    if (in[0] == in[1]) {
        return T(1.);
    }
    return in[0] / in[1];
}

std::string print_my_pdiv(const std::vector<std::string> &in)
{
    return "(" + in[0] + "/" + in[1] + ")";
}

/*--------------------------------------------------------------------------
 *                                  UNARY FUNCTIONS
 *------------------------------------------------------------------------**/
// sine
template <typename T, f_enabler<T> = 0>
T my_sin(const std::vector<T> &in)
{
    return sin(in[0]);
}

std::string print_my_sin(const std::vector<std::string> &in)
{
    return "sin(" + in[0] + ")";
}

// cosine
template <typename T, f_enabler<T> = 0>
T my_cos(const std::vector<T> &in)
{
    return cos(in[0]);
}

std::string print_my_cos(const std::vector<std::string> &in)
{
    return "cos(" + in[0] + ")";
}

// logarithm
template <typename T, f_enabler<T> = 0>
T my_log(const std::vector<T> &in)
{
    return audi::log(in[0]);
}

std::string print_my_log(const std::vector<std::string> &in)
{
    return "log(" + in[0] + ")";
}

// exponential
template <typename T, f_enabler<T> = 0>
T my_exp(const std::vector<T> &in)
{
    return audi::exp(in[0]);
}

std::string print_my_exp(const std::vector<std::string> &in)
{
    return "exp(" + in[0] + ")";
}

} // namespace dcgp

#endif // DCGP_WRAPPED_FUNCTIONS_H
