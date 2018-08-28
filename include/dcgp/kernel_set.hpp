#ifndef DCGP_kernel_set_H
#define DCGP_kernel_set_H

#include <audi/gdual.hpp>
#include <vector>

#include <dcgp/kernel.hpp>
#include <dcgp/wrapped_functions.hpp>

namespace dcgp
{

/// Function set
/**
 * This class is provided as an helper to construct the std::vector<kernel<T>>
 * that is requested to form a dcgp::expression<T>. Once constructed, a call to its
 * parenthesis operator will return the std::vector containing the requested kernels
 *
 * @tparam T The type of the functions output (and inputs)
 */
template <typename T>
class kernel_set
{
public:
    /// Constructor
    /**
     * Default constructor
     */
    kernel_set() : m_kernels(){};

    /// Constructor
    /**
     * Constructs a kernel set that can be used in a dCGP expression
     *
     * @param[in] list an std::vector of strings containing the function names (e.g., "sum")
     */
    kernel_set(const std::vector<std::string> &list)
    {
        for (auto kernel_name : list) {
            push_back(kernel_name);
        }
    }

    /// Adds a kernel to the set
    /**
     * Constructs a kernel<T> given a string containing the function name, and
     * inserts it into the std::vector
     *
     * @param[in] kernel_name a string containing the function name
     *
     * @throw std::invalid_argument if the function is not implemented
     */
    void push_back(std::string kernel_name)
    {
        if (kernel_name == "sum")
            m_kernels.emplace_back(my_sum<T>, print_my_sum, kernel_name);
        else if (kernel_name == "diff")
            m_kernels.emplace_back(my_diff<T>, print_my_diff, kernel_name);
        else if (kernel_name == "mul")
            m_kernels.emplace_back(my_mul<T>, print_my_mul, kernel_name);
        else if (kernel_name == "div")
            m_kernels.emplace_back(my_div<T>, print_my_div, kernel_name);
        else if (kernel_name == "pdiv")
            m_kernels.emplace_back(my_pdiv<T>, print_my_pdiv, kernel_name);
        else if (kernel_name == "sig")
            m_kernels.emplace_back(my_sig<T>, print_my_sig, kernel_name);
        else if (kernel_name == "tanh")
            m_kernels.emplace_back(my_tanh<T>, print_my_tanh, kernel_name);
        else if (kernel_name == "ReLu")
            m_kernels.emplace_back(my_relu<T>, print_my_relu, kernel_name);
        else if (kernel_name == "sin")
            m_kernels.emplace_back(my_sin<T>, print_my_sin, kernel_name);
        else if (kernel_name == "cos")
            m_kernels.emplace_back(my_cos<T>, print_my_cos, kernel_name);
        else if (kernel_name == "log")
            m_kernels.emplace_back(my_log<T>, print_my_log, kernel_name);
        else if (kernel_name == "exp")
            m_kernels.emplace_back(my_exp<T>, print_my_exp, kernel_name);
        else
            throw std::invalid_argument("Unimplemented function " + kernel_name);
    }

    /// Adds a kernel to the set
    /**
     * Inserts the given kernel<T> into the std::vector
     *
     * @param[in] kernel the dcgp::kernel<T> to add
     */
    void push_back(const dcgp::kernel<T> &kernel)
    {
        m_kernels.push_back(kernel);
    }

    /// Clears the kernel set
    /**
     * Removes all the elements from the std::vector containing the kernels
     */
    void clear()
    {
        m_kernels.clear();
    }

    /// Overloaded function call operator
    /**
     * Returns the std::vector containing the kernels
     */
    std::vector<dcgp::kernel<T>> operator()() const
    {
        return m_kernels;
    }

    /// Overloaded stream operator
    /**
     * Will stream the function names
     *
     * @param[in,out] os target stream
     * @param[in] d kernel_set argument
     *
     * @return reference to \p os
     */
    friend std::ostream &operator<<(std::ostream &os, const kernel_set<T> &d)
    {
        stream(os, d());
        return os;
    }

    // Overloaded subscript operator
    /**
     * Returns the kernel at a specified index in the std::vector
     *
     * @param[in] idx index of the required kernel
     *
     * @return the kernel at the specified index
     */
    dcgp::kernel<T> operator[](const typename std::vector<dcgp::kernel<T>>::size_type idx) const
    {
        return m_kernels[idx];
    }

private:
    // vector of functions
    std::vector<dcgp::kernel<T>> m_kernels;
};

} // end of namespace dcgp

#endif // DCGP_kernel_set_H
