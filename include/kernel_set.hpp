#ifndef DCGP_kernel_set_H
#define DCGP_kernel_set_H

#include <vector>
#include <audi/gdual.hpp>

#include "kernel.hpp"
#include "wrapped_functions.hpp"

namespace dcgp {

/// Function set
/**
 * This class is provided as an helper to construct the std::vector<kernel<T>>
 * that is requested to form a dcgp::expression<T>. Once constructed, a call to its
 * parenthesis operator will return the std::vector containing the requested kernels
 */
template<typename T>
class kernel_set
{
public:
    kernel_set() : m_kernels() {};
    kernel_set(const std::vector<std::string>& list)
    {
        for (auto kernel_name : list)
        {
            push_back(kernel_name);
        }
    };

    void push_back(std::string kernel_name)
    {
        if (kernel_name=="sum")
            m_kernels.emplace_back(my_sum<T>, print_my_sum, kernel_name);
        else if (kernel_name=="diff")
            m_kernels.emplace_back(my_diff<T>, print_my_diff, kernel_name);
        else if (kernel_name=="mul")
            m_kernels.emplace_back(my_mul<T>,print_my_mul, kernel_name);
        else if (kernel_name=="div")
            m_kernels.emplace_back(my_div<T>,print_my_div, kernel_name);
        else if (kernel_name=="sig")
            m_kernels.emplace_back(my_sig<T>,print_my_sig, kernel_name);
        else if (kernel_name=="sin")
            m_kernels.emplace_back(my_sin<T>,print_my_sin, kernel_name);
        else if (kernel_name=="cos")
            m_kernels.emplace_back(my_cos<T>,print_my_cos, kernel_name);
        else if (kernel_name=="log")
            m_kernels.emplace_back(my_log<T>,print_my_log, kernel_name);
        else if (kernel_name=="exp")
            m_kernels.emplace_back(my_exp<T>,print_my_exp, kernel_name);
        else
            throw std::invalid_argument("Unimplemented function " + kernel_name);
    };

    void push_back(const dcgp::kernel<T>& kernel)
    {
        m_kernels.push_back(kernel);
    };

    void clear()
    {
        m_kernels.clear();
    };

    std::vector<dcgp::kernel<T>> operator()() const
    {
        return m_kernels;
    };

    friend std::ostream& operator<<(std::ostream& os, const kernel_set<T>& d)
    {
        stream(os, d());
        return os;
    }

    dcgp::kernel<T> operator[] (const typename std::vector<dcgp::kernel<T>>::size_type idx) const {
        return m_kernels[idx];
    }

private:
    std::vector<dcgp::kernel<T>> m_kernels;
};

} // end of namespace dcgp

#endif // DCGP_kernel_set_H
