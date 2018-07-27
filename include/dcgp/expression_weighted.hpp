#ifndef DCGP_EXPRESSION_WEIGHTED_H
#define DCGP_EXPRESSION_WEIGHTED_H

#include <audi/audi.hpp>
#include <initializer_list>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/kernel.hpp>
#include <dcgp/type_traits.hpp>

namespace dcgp
{

/// A weighted dCGP expression
/**
 * This class represents a mathematical expression as encoded using CGP with the
 * addition of weights on the connections. It contains algorithms to compute the
 * value (numerical and symbolical) of the expression and its derivatives, as well
 * as to mutate the expression.
 *
 * @tparam T expression type. Can be double, or a gdual type.
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
template <typename T>
class expression_weighted : public expression<T>
{

private:
    // SFINAE dust
    template <typename U>
    using functor_enabler = typename std::enable_if<
        std::is_same<U, double>::value || is_gdual<T>::value || std::is_same<U, std::string>::value, int>::type;

public:
    /// Constructor
    /** Constructs a dCGP expression
     *
     * @param[in] n number of inputs (independent variables)
     * @param[in] m number of outputs (dependent variables)
     * @param[in] r number of rows of the dCGP
     * @param[in] c number of columns of the dCGP
     * @param[in] l number of levels-back allowed for the dCGP
     * @param[in] arity arity of the basis functions
     * @param[in] f function set. An std::vector of dcgp::kernel<expression::type>
     * @param[in] seed seed for the random number generator (initial expression and mutations depend on this)
     */
    expression_weighted(unsigned int n,           // n. inputs
                        unsigned int m,           // n. outputs
                        unsigned int r,           // n. rows
                        unsigned int c,           // n. columns
                        unsigned int l,           // n. levels-back
                        unsigned int arity,       // basis functions' arity
                        std::vector<kernel<T>> f, // functions
                        unsigned int seed         // seed for the pseudo-random numbers
                        )
        : expression<T>(n, m, r, c, l, arity, f, seed), m_weights(r * c * arity, T(1.))
    {
        for (auto i = 0u; i < r * c; ++i) {
            for (auto j = 0u; j < arity; ++j) {
                m_weights_symbols.push_back("w" + std::to_string(i + n) + "_" + std::to_string(j));
            }
        }
    }

    /// Evaluates the dCGP expression
    /**
     * This evaluates the dCGP expression. According to the template parameter
     * it will compute the value (double) the Taylor expansion (gdual) or a symbolic
     * representation (std::string). Any other type will result in a compilation-time
     * error (SFINAE).
     *
     * @param[in] in an std::vector containing the values where the dCGP expression has
     * to be computed (doubles, gduals or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename U, functor_enabler<U> = 0>
    std::vector<U> operator()(const std::vector<U> &in) const
    {
        if (in.size() != this->get_n()) {
            throw std::invalid_argument("Input size is incompatible");
        }
        std::vector<U> retval(this->get_m());
        std::unordered_map<unsigned int, U> node;
        std::vector<U> function_in(this->get_arity());
        for (auto i : this->get_active_nodes()) {
            if (i < this->get_n()) {
                node[i] = in[i];
            } else {
                unsigned int idx
                    = (i - this->get_n()) * (this->get_arity() + 1); // position in the chromosome of the current node
                unsigned int weight_idx = (i - this->get_n()) * this->get_arity();
                for (auto j = 0u; j < this->get_arity(); ++j) {
                    function_in[j] = node[this->get()[idx + j + 1]];
                }
                node[i] = kernel_call(function_in, idx, weight_idx);
            }
        }
        for (auto i = 0u; i < this->get_m(); ++i) {
            retval[i] = node[this->get()[(this->get_rows() * this->get_cols()) * (this->get_arity() + 1) + i]];
        }
        return retval;
    }

    /// Evaluates the dCGP expression
    /**
     * This evaluates the dCGP expression. According to the template parameter
     * it will compute the value (double) the Taylor expansion (gdual) or a symbolic
     * representation (std::string). Any other type will result in a compilation-time
     * error (SFINAE). This is identical to the other overload and is provided only
     * for convenience
     *
     * @param[in] in an initializer list containing the values where the dCGP expression has
     * to be computed (doubles, gduals or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename U, functor_enabler<U> = 0>
    std::vector<U> operator()(const std::initializer_list<U> &in) const
    {
        std::vector<U> dummy(in);
        return (*this)(dummy);
    }

    /// Overloaded stream operator
    /**
     * Will return a formatted string containing a human readable representation
     * of the class
     *
     * @return std::string containing a human-readable representation of the problem.
     */
    friend std::ostream &operator<<(std::ostream &os, const expression_weighted &d)
    {
        audi::stream(os, "d-CGP Expression:\n");
        audi::stream(os, "\tNumber of inputs:\t\t", d.get_n(), '\n');
        audi::stream(os, "\tNumber of outputs:\t\t", d.get_m(), '\n');
        audi::stream(os, "\tNumber of rows:\t\t\t", d.get_rows(), '\n');
        audi::stream(os, "\tNumber of columns:\t\t", d.get_cols(), '\n');
        audi::stream(os, "\tNumber of levels-back allowed:\t", d.get_levels_back(), '\n');
        audi::stream(os, "\tBasis function arity:\t\t", d.get_arity(), '\n');
        audi::stream(os, "\n\tResulting lower bounds:\t", d.get_lb());
        audi::stream(os, "\n\tResulting upper bounds:\t", d.get_ub(), '\n');
        audi::stream(os, "\n\tCurrent expression (encoded):\t", d.get(), '\n');
        audi::stream(os, "\tActive nodes:\t\t\t", d.get_active_nodes(), '\n');
        audi::stream(os, "\tActive genes:\t\t\t", d.get_active_genes(), '\n');
        audi::stream(os, "\n\tFunction set:\t\t\t", d.get_f(), '\n');
        audi::stream(os, "\n\tWeights:\t\t\t", d.m_weights, '\n');

        return os;
    }

    /// Sets a weight
    /**
     * Sets a connection weight to a new value
     *
     * @param[in] node_id the id of the node whose weight is being set (convention adopted for node numbering
     * http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf)
     * @param[in] input_id the id of the node input (0 for the first one up to arity-1)
     * @param[in] w the new value of the weight
     *
     * @throws std::invalid_argument if the node_id or input_id are not valid
     */
    void set_weight(typename std::vector<T>::size_type node_id, typename std::vector<T>::size_type input_id, const T &w)
    {
        if (node_id < this->get_n() || node_id >= this->get_n() + this->get_rows() * this->get_cols()) {
            throw std::invalid_argument("Requested node id does not exist");
        }
        if (input_id >= this->get_arity()) {
            throw std::invalid_argument("Requested input exceeds the function arity");
        }
        auto idx = (node_id - this->get_n()) * this->get_arity() + input_id;
        m_weights[idx] = w;
    }

    /// Gets a weight
    /**
     * Gets the value of a connection weight
     *
     * @param[in] node_id the id of the node (convention adopted for node numbering
     * http://ppsn2014.ijs.si/files/slides/ppsn2014-tutorial3-miller.pdf)
     * @param[in] input_id the id of the node input (0 for the first one up to arity-1)
     *
     * @return the value of the weight
     *
     * @throws std::invalid_argument if the node_id or input_id are not valid
     */
    T get_weight(typename std::vector<T>::size_type node_id, typename std::vector<T>::size_type input_id)
    {
        if (node_id < this->get_n() || node_id >= this->get_n() + this->get_rows() * this->get_cols()) {
            throw std::invalid_argument("Requested node id does not exist");
        }
        if (input_id >= this->get_arity()) {
            throw std::invalid_argument("Requested input exceeds the function arity");
        }
        auto idx = (node_id - this->get_n()) * this->get_arity() + input_id;
        return m_weights[idx];
    }

    /// Sets all weights
    /**
     * Sets all the connection weights at once
     *
     * @param[in] ws an std::vector containing all the weights to set
     *
     * @throws std::invalid_argument if the input vector dimension is not valid (r*c*arity)
     */
    void set_weights(const std::vector<T> &ws)
    {
        if (ws.size() != m_weights.size()) {
            throw std::invalid_argument("The vector of weights has the wrong dimension");
        }
        m_weights = ws;
    }

    /// Gets the weights
    /**
     * Gets the values of all the weights
     *
     * @return an std::vector containing all the weights
     */
    const std::vector<T> &get_weights() const
    {
        return m_weights;
    }

private:

    // For numeric computations
    template <typename U, typename std::enable_if<std::is_same<U, double>::value || is_gdual<U>::value, int>::type = 0>
    U kernel_call(std::vector<U> &function_in, unsigned int idx, unsigned int weight_idx) const
    {
        for (auto j = 0u; j < this->get_arity(); ++j) {
            function_in[j] = function_in[j] * m_weights[weight_idx + j];
        }
        return this->get_f()[this->get()[idx]](function_in);
    }

    // For the symbolic expression
    template <typename U, typename std::enable_if<std::is_same<U, std::string>::value, int>::type = 0>
    U kernel_call(std::vector<U> &function_in, unsigned int idx, unsigned int weight_idx) const
    {
        for (auto j = 0u; j < this->get_arity(); ++j) {
            function_in[j] = "(" + m_weights_symbols[weight_idx + j] + "*" + function_in[j] + ")";
        }
        return this->get_f()[this->get()[idx]](function_in);
    }

    std::vector<T> m_weights;
    std::vector<std::string> m_weights_symbols;
};

} // end of namespace dcgp

#endif // DCGP_EXPRESSION_H
