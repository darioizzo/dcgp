#include <audi/gdual.hpp>

#include "function_set.h"
#include "wrapped_functions.h"

namespace dcgp {

function_set::function_set() : m_functions() {};

function_set::function_set(const std::vector<std::string>& list)
{
    for (auto function_name : list)
    {
        push_back(function_name);
    }
}

void function_set::push_back(const std::string& function_name)
{
    if (function_name=="sum")
        m_functions.emplace_back(my_sum<double>, my_sum<audi::gdual>, print_my_sum, function_name);
    else if (function_name=="diff")
        m_functions.emplace_back(my_diff<double>,my_diff<audi::gdual>,print_my_diff, function_name);
    else if (function_name=="mul")
        m_functions.emplace_back(my_mul<double>,my_mul<audi::gdual>,print_my_mul, function_name);
    else if (function_name=="div")
        m_functions.emplace_back(my_div<double>,my_div<audi::gdual>,print_my_div, function_name);
    else 
        throw std::invalid_argument("Unimplemented function " + function_name);
}

void function_set::clear()
{
    m_functions.clear();
}

std::vector<basis_function> function_set::operator()() const
{
    return m_functions;
}


} // end of namespace dcgp
