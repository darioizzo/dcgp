#ifndef DCGP_FUNCTION_SET_H
#define DCGP_FUNCTION_SET_H

#include <vector>
#include <audi/gdual.hpp>

#include "basis_function.hpp"
#include "wrapped_functions.hpp"

namespace dcgp {

/// Function set
/**
 * Contains,
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class function_set
{
public:
    function_set() : m_functions() {};
    function_set(const std::vector<std::string>& list)
    {
        for (auto function_name : list)
        {
            push_back(function_name);
        }
    };

    void push_back(const std::string& function_name)
    {
        if (function_name=="sum")
            m_functions.emplace_back(my_sum<double>, my_sum<audi::gdual>, print_my_sum, function_name);
        else if (function_name=="diff")
            m_functions.emplace_back(my_diff<double>,my_diff<audi::gdual>,print_my_diff, function_name);
        else if (function_name=="mul")
            m_functions.emplace_back(my_mul<double>,my_mul<audi::gdual>,print_my_mul, function_name);
        else if (function_name=="div")
            m_functions.emplace_back(my_div<double>,my_div<audi::gdual>,print_my_div, function_name);
        else if (function_name=="sig")
            m_functions.emplace_back(my_sig<double>,my_sig<audi::gdual>,print_my_sig, function_name);
        else if (function_name=="sin")
            m_functions.emplace_back(my_sin<double>,my_sin<audi::gdual>,print_my_sin, function_name);
        else if (function_name=="log")
            m_functions.emplace_back(my_log<double>,my_log<audi::gdual>,print_my_log, function_name);
        else if (function_name=="exp")
            m_functions.emplace_back(my_exp<double>,my_exp<audi::gdual>,print_my_exp, function_name);
        else
            throw std::invalid_argument("Unimplemented function " + function_name);
    };

    void clear()
    {
        m_functions.clear();
    };

    std::vector<dcgp::basis_function> operator()() const
    {
        return m_functions;
    };
private:
    std::vector<dcgp::basis_function> m_functions;
};

} // end of namespace dcgp

#endif // DCGP_FUNCTION_SET_H
