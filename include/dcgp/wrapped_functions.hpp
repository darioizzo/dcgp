#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <audi/audi.hpp>
#include <audi/functions.hpp>
#include <cmath>
#include <string>
#include <type_traits>
#include <vector>

#include <dcgp/config.hpp>
#include <dcgp/function.hpp>
#include <dcgp/type_traits.hpp>

#define DCGP_S11N_FUNCTION_EXPORT_KEY_MULTI(f)                                                                         \
    DCGP_S11N_FUNCTION_EXPORT_KEY(dcgp_##f##_double, dcgp::f##_func<double>, double, const std::vector<double> &)      \
    DCGP_S11N_FUNCTION_EXPORT_KEY(dcgp_##f##_gdual_d, dcgp::f##_func<audi::gdual_d>, audi::gdual_d,                    \
                                  const std::vector<audi::gdual_d> &)                                                  \
    DCGP_S11N_FUNCTION_EXPORT_KEY(dcgp_##f##_gdual_v, dcgp::f##_func<audi::gdual_v>, audi::gdual_v,                    \
                                  const std::vector<audi::gdual_v> &)

#define DCGP_S11N_FUNCTION_IMPLEMENT_MULTI(f)                                                                          \
    DCGP_S11N_FUNCTION_IMPLEMENT(dcgp_##f##_double, dcgp::f##_func<double>, double, const std::vector<double> &)       \
    DCGP_S11N_FUNCTION_IMPLEMENT(dcgp_##f##_gdual_d, dcgp::f##_func<audi::gdual_d>, audi::gdual_d,                     \
                                 const std::vector<audi::gdual_d> &)                                                   \
    DCGP_S11N_FUNCTION_IMPLEMENT(dcgp_##f##_gdual_v, dcgp::f##_func<audi::gdual_v>, audi::gdual_v,                     \
                                 const std::vector<audi::gdual_v> &)

namespace dcgp
{

// SFINAE dust (to hide under the carpet). Its used to enable the templated
// version of the various functions that can construct a kernel object. Only for
// double and a gdual type Complex could also be allowed.
template <typename T>
using f_enabler = typename std::enable_if<std::is_same<T, double>::value || is_gdual<T>::value, int>::type;

/*--------------------------------------------------------------------------
 *                              N-ARITY FUNCTIONS
 *------------------------------------------------------------------------**/

template <typename T, f_enabler<T> = 0>
struct my_diff_func {
    T operator()(const std::vector<T> &in) const
    {
        T retval(in[0]);
        for (auto i = 1u; i < in.size(); ++i) {
            retval -= in[i];
        }
        return retval;
    }
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

template <typename T>
inline constexpr auto my_diff = my_diff_func<T>{};

inline std::string print_my_diff(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "-" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
inline T my_mul(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval *= in[i];
    }
    return retval;
}

inline std::string print_my_mul(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "*" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
inline T my_div(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval /= in[i];
    }
    return retval;
}

inline std::string print_my_div(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "/" + in[i];
    }
    return "(" + retval + ")";
}

// Protected divide function (double overload):
template <typename T, typename std::enable_if<std::is_same<T, double>::value, int>::type = 0>
inline T my_pdiv(const std::vector<T> &in)
{
    T retval(in[0]);
    T tmpval(in[1]);

    for (auto i = 2u; i < in.size(); ++i) {
        tmpval *= in[i];
    }

    retval /= tmpval;

    if (std::isfinite(retval)) {
        return retval;
    }

    return 1.;
}

// Protected divide function (gdual overload):
// this will throw a compiler error when used.
// The pdiv is only available as a double type, for use in CGP.
// Because the gradients created when using gdual are mathematically invalid.
template <typename T, typename std::enable_if<is_gdual<T>::value, int>::type = 0>
inline T my_pdiv(const std::vector<T> &)
{
    throw std::invalid_argument("The protected division is not supported for gdual types.");
}

inline std::string print_my_pdiv(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "/" + in[i];
    }

    return "(" + retval + ")";
}

/*--------------------------------------------------------------------------
 *                            Suitable for dCGPANN
 *------------------------------------------------------------------------**/

// sigmoid function: 1 / (1 + exp(- (a + b + c + d+ .. + ))
template <typename T, f_enabler<T> = 0>
inline T my_sig(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return 1. / (1. + audi::exp(-retval));
}

inline std::string print_my_sig(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "sig(" + retval + ")";
}

// tanh function:
template <typename T, f_enabler<T> = 0>
inline T my_tanh(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return audi::tanh(retval);
}

inline std::string print_my_tanh(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "tanh(" + retval + ")";
}

// ReLu function (double overload):
template <typename T, typename std::enable_if<std::is_same<T, double>::value, int>::type = 0>
inline T my_relu(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    if (retval < 0) {
        retval = T(0.);
    }
    return retval;
}

// ReLu function (gdual overload):
template <typename T, typename std::enable_if<is_gdual<T>::value, int>::type = 0>
inline T my_relu(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    if (retval.constant_cf() < T(0.).constant_cf()) {
        retval = T(0.);
    }
    return retval;
}

inline std::string print_my_relu(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "ReLu(" + retval + ")";
}

// Exponential linear unit (ELU) function (double overload):
template <typename T, typename std::enable_if<std::is_same<T, double>::value, int>::type = 0>
inline T my_elu(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    if (retval < 0) {
        retval = audi::exp(retval) - T(1.);
    }
    return retval;
}

// Exponential linear unit (ELU) function (gdual overload):
template <typename T, typename std::enable_if<is_gdual<T>::value, int>::type = 0>
inline T my_elu(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    if (retval.constant_cf() < T(0.).constant_cf()) {
        retval = audi::exp(retval) - T(1.);
    }
    return retval;
}

inline std::string print_my_elu(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "ELU(" + retval + ")";
}

// Inverse square root function: x / sqrt(1+x^2):
template <typename T, f_enabler<T> = 0>
inline T my_isru(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return retval / (audi::sqrt(1 + retval * retval));
}

inline std::string print_my_isru(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "ISRU(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
inline T my_sum(const std::vector<T> &in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += in[i];
    }
    return retval;
}

inline std::string print_my_sum(const std::vector<std::string> &in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval += "+" + in[i];
    }
    return "(" + retval + ")";
}

/*--------------------------------------------------------------------------
 *                               UNARY FUNCTIONS
 *------------------------------------------------------------------------**/
// sine
template <typename T, f_enabler<T> = 0>
inline T my_sin(const std::vector<T> &in)
{
    return sin(in[0]);
}

inline std::string print_my_sin(const std::vector<std::string> &in)
{
    return "sin(" + in[0] + ")";
}

// cosine
template <typename T, f_enabler<T> = 0>
inline T my_cos(const std::vector<T> &in)
{
    return cos(in[0]);
}

inline std::string print_my_cos(const std::vector<std::string> &in)
{
    return "cos(" + in[0] + ")";
}

// logarithm
template <typename T, f_enabler<T> = 0>
inline T my_log(const std::vector<T> &in)
{
    return audi::log(in[0]);
}

inline std::string print_my_log(const std::vector<std::string> &in)
{
    return "log(" + in[0] + ")";
}

// exponential (unary)
// This exponential discards all inputs except the first one
template <typename T, f_enabler<T> = 0>
inline T my_exp(const std::vector<T> &in)
{
    return audi::exp(in[0]);
}

inline std::string print_my_exp(const std::vector<std::string> &in)
{
    return "exp(" + in[0] + ")";
}

// gaussian (unary)
// This gaussian discards all inputs except the first one
template <typename T, f_enabler<T> = 0>
inline T my_gaussian(const std::vector<T> &in)
{
    return audi::exp(-in[0] * in[0]);
}

inline std::string print_my_gaussian(const std::vector<std::string> &in)
{
    return "exp(-" + in[0] + "**2)";
}

// sqrt (unary)
// This square root discards all inputs except the first one
template <typename T, f_enabler<T> = 0>
struct my_sqrt_func {
    T operator()(const std::vector<T> &in) const
    {
        return audi::sqrt(in[0]);
    }
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

template <typename T>
inline constexpr auto my_sqrt = my_sqrt_func<T>{};

inline std::string print_my_sqrt(const std::vector<std::string> &in)
{
    return "sqrt(" + in[0] + ")";
}

} // namespace dcgp

DCGP_S11N_FUNCTION_EXPORT_KEY_MULTI(my_diff)
DCGP_S11N_FUNCTION_EXPORT_KEY_MULTI(my_sqrt)

#endif // DCGP_WRAPPED_FUNCTIONS_H
