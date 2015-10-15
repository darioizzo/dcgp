#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <string>
#include <vector>
#include <audi/gdual.hpp>
#include <audi/functions.hpp>

namespace dcgp {

// SFINAE dust (to hide under the carpet)
template <typename T>
using f_enabler = typename std::enable_if<std::is_same<T,double>::value || std::is_same<T,audi::gdual>::value, int>::type;

// Allows to overload in templates std functions with audi functions
using namespace std;
using namespace audi;

/*--------------------------------------------------------------------------
*                                  BINARY FUNCTIONS
*------------------------------------------------------------------------**/
template <typename T, f_enabler<T> = 0>
T my_sum(const std::vector<T>& in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+=in[i];
    }
    return retval;
}

std::string print_my_sum(const std::vector<std::string>& in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+= "+" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
T my_diff(const std::vector<T>& in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval-=in[i];
    }
    return retval;
}

std::string print_my_diff(const std::vector<std::string>& in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+= "-" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
T my_mul(const std::vector<T>& in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval*=in[i];
    }
    return retval;
}

std::string print_my_mul(const std::vector<std::string>& in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+= "*" + in[i];
    }
    return "(" + retval + ")";
}

template <typename T, f_enabler<T> = 0>
T my_div(const std::vector<T>& in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval/=in[i];
    }
    return retval;
}

std::string print_my_div(const std::vector<std::string>& in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+= "/" + in[i];
    }
    return "(" + retval + ")";
}

// sigmoid function: 1 / (1 + exp(- (a + b + c + d+ .. + ))
template <typename T, f_enabler<T> = 0>
T my_sig(const std::vector<T>& in)
{
    T retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+=in[i];
    }
    return 1. / (1. + exp(-retval));
}

std::string print_my_sig(const std::vector<std::string>& in)
{
    std::string retval(in[0]);
    for (auto i = 1u; i < in.size(); ++i) {
        retval+= "+" + in[i];
    }
    return "sigmoid(" + retval + ")";
}

/*--------------------------------------------------------------------------
*                                  UNARY FUNCTIONS
*------------------------------------------------------------------------**/
// sine
template <typename T, f_enabler<T> = 0>
T my_sin(const std::vector<T>& in)
{
    return sin(in[0]);
}

std::string print_my_sin(const std::vector<std::string>& in)
{
    return "sin(" + in[0] + ")";
}

// logarithm
template <typename T, f_enabler<T> = 0>
T my_log(const std::vector<T>& in)
{
    return log(in[0]);
}

std::string print_my_log(const std::vector<std::string>& in)
{
    return "log(" + in[0] + ")";
}

// logarithm
template <typename T, f_enabler<T> = 0>
T my_exp(const std::vector<T>& in)
{
    return exp(in[0]);
}

std::string print_my_exp(const std::vector<std::string>& in)
{
    return "exp(" + in[0] + ")";
}

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H