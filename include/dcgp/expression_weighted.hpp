#ifndef DCGP_EXPRESSION_WEIGHTED_H
#define DCGP_EXPRESSION_WEIGHTED_H

#include <initializer_list>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <audi/audi.hpp>

#include <dcgp/config.hpp>
#include <dcgp/expression.hpp>
#include <dcgp/kernel.hpp>
#include <dcgp/s11n.hpp>
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
    /** Constructs a weighted dCGP expression.
     *
     * @param[in] n number of inputs (independent variables).
     * @param[in] m number of outputs (dependent variables).
     * @param[in] r number of rows of the weighted dCGP.
     * @param[in] c number of columns of the weighted dCGP.
     * @param[in] l number of levels-back allowed for the weighted dCGP.
     * @param[in] arity arities of the basis functions for each column.
     * @param[in] f function set. An std::vector of dcgp::kernel<expression::type>.
     * @param[in] seed seed for the random number generator (initial expression and mutations depend on this).
     */
    expression_weighted(unsigned n = 1u,                                     // n. inputs
                        unsigned m = 1u,                                     // n. outputs
                        unsigned r = 1u,                                     // n. rows
                        unsigned c = 1u,                                     // n. columns
                        unsigned l = 1u,                                     // n. levels-back
                        unsigned arity = 2u,                                 // basis functions' arity
                        std::vector<kernel<T>> f = kernel_set<T>({"sum"})(), // functions
                        unsigned seed = dcgp::random_device::next())
        : expression<T>(n, m, r, c, l, std::vector<unsigned>(c, arity), f, 0u, seed)
    {
        // Default initialization of weights to 1.
        unsigned n_connections = std::accumulate(this->get_arity().begin(), this->get_arity().end(), 0u) * r;
        m_weights = std::vector<T>(n_connections, T(1.));

        // Filling in the symbols for the weights and biases
        for (auto node_id = n; node_id < r * c + n; ++node_id) {
            for (auto j = 0u; j < this->_get_arity(node_id); ++j) {
                m_weights_symbols.push_back("w" + std::to_string(node_id) + "_" + std::to_string(j));
            }
        }
    }

    /// Constructor
    /** Constructs a weighted dCGP expression
     *
     * @param[in] n number of inputs (independent variables).
     * @param[in] m number of outputs (dependent variables).
     * @param[in] r number of rows of the weighted dCGP.
     * @param[in] c number of columns of the weighted dCGP.
     * @param[in] l number of levels-back allowed for the weighted dCGP.
     * @param[in] arities arities of the basis functions (per each column).
     * @param[in] f function set. An std::vector of dcgp::kernel<expression::type>.
     * @param[in] seed seed for the random number generator (initial expression and mutations depend on this).
     */
    expression_weighted(unsigned n,                    // n. inputs
                        unsigned m,                    // n. outputs
                        unsigned r,                    // n. rows
                        unsigned c,                    // n. columns
                        unsigned l,                    // n. levels-back
                        std::vector<unsigned> arities, // basis functions' arity
                        std::vector<kernel<T>> f,      // functions
                        unsigned seed = dcgp::random_device::next())
        : expression<T>(n, m, r, c, l, arities, f, 0u, seed)
    {
        // Default initialization of weights to 1.
        unsigned n_connections = std::accumulate(this->get_arity().begin(), this->get_arity().end(), 0u) * r;
        m_weights = std::vector<T>(n_connections, T(1.));

        // Filling in the symbols for the weights and biases.
        for (auto node_id = n; node_id < r * c + n; ++node_id) {
            for (auto j = 0u; j < this->_get_arity(node_id); ++j) {
                m_weights_symbols.push_back("w" + std::to_string(node_id) + "_" + std::to_string(j));
            }
        }
    }

    /// Evaluates the dCGP-weighted expression
    /**
     * This evaluates the dCGP-weighted expression. This method overrides the base class
     * method. NOTE we cannot template this and the following function as they are virtual in the base class.
     *
     * @param[in] in std::vector containing the values where the dCGP-weighted expression has
     * to be computed
     *
     * @return The value of the output (an std::vector)
     */
    std::vector<T> operator()(const std::vector<T> &in) const
    {
        if (in.size() != this->get_n()) {
            throw std::invalid_argument("Input size is incompatible");
        }
        std::vector<T> retval(this->get_m());
        std::vector<T> node(this->get_n() + this->get_r() * this->get_c());
        std::vector<T> function_in;
        for (auto node_id : this->get_active_nodes()) {
            if (node_id < this->get_n()) {
                node[node_id] = in[node_id];
            } else {
                unsigned arity = this->_get_arity(node_id);
                function_in.resize(arity);
                // position in the chromosome of the current node
                unsigned g_idx = this->get_gene_idx()[node_id];
                // starting position in m_weights of the weights relative to the node
                unsigned w_idx = g_idx - (node_id - this->get_n());
                for (unsigned j = 0u; j < this->_get_arity(node_id); ++j) {
                    function_in[j] = node[this->get()[g_idx + j + 1]];
                }
                node[node_id] = kernel_call(function_in, g_idx, node_id, w_idx);
            }
        }
        for (auto i = 0u; i < this->get_m(); ++i) {
            retval[i] = node[this->get()[this->get().size() - this->get_m() + i]];
        }
        return retval;
    }

    /// Evaluates the dCGP-weighted expression
    /**
     * This evaluates the dCGP-weighted expression. This method overrides the base class
     * method.
     *
     * @param[in] in an std::vector containing the symbol names.
     *
     * @return The symbolic value of the output (an std::vector)
     */
    std::vector<std::string> operator()(const std::vector<std::string> &in) const
    {
        if (in.size() != this->get_n()) {
            throw std::invalid_argument("Input size is incompatible");
        }
        std::vector<std::string> retval(this->get_m());
        std::vector<std::string> node(this->get_n() + this->get_r() * this->get_c());
        std::vector<std::string> function_in;
        for (auto node_id : this->get_active_nodes()) {
            if (node_id < this->get_n()) {
                node[node_id] = in[node_id];
            } else {
                unsigned arity = this->_get_arity(node_id);
                function_in.resize(arity);
                // position in the chromosome of the current node
                unsigned g_idx = this->get_gene_idx()[node_id];
                // starting position in m_weights of the weights relative to the node
                unsigned w_idx = g_idx - (node_id - this->get_n());
                for (unsigned j = 0u; j < this->_get_arity(node_id); ++j) {
                    function_in[j] = node[this->get()[g_idx + j + 1]];
                }
                node[node_id] = kernel_call(function_in, g_idx, node_id, w_idx);
            }
        }
        for (auto i = 0u; i < this->get_m(); ++i) {
            retval[i] = node[this->get()[this->get().size() - this->get_m() + i]];
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
        audi::stream(os, "\tNumber of rows:\t\t\t", d.get_r(), '\n');
        audi::stream(os, "\tNumber of columns:\t\t", d.get_c(), '\n');
        audi::stream(os, "\tNumber of levels-back allowed:\t", d.get_l(), '\n');
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
    void set_weight(unsigned node_id, unsigned input_id, const T &w)
    {
        if (node_id < this->get_n() || node_id >= this->get_n() + this->get_r() * this->get_c()) {
            throw std::invalid_argument("Requested node id does not exist");
        }
        if (input_id >= this->_get_arity(node_id)) {
            throw std::invalid_argument("Requested input exceeds the function arity");
        }
        // index of the node in the weight vector
        auto idx = this->get_gene_idx()[node_id] - (node_id - this->get_n());
        m_weights[idx + input_id] = w;
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
    T get_weight(unsigned node_id, unsigned input_id) const
    {
        if (node_id < this->get_n() || node_id >= this->get_n() + this->get_r() * this->get_c()) {
            throw std::invalid_argument(
                "Requested node id does not exist or does not have a weight (e.g. input nodes)");
        }
        if (input_id >= this->_get_arity(node_id)) {
            throw std::invalid_argument("Requested input exceeds the function arity");
        }

        auto idx = this->get_gene_idx()[node_id] - (node_id - this->get_n());
        return m_weights[idx + input_id];
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

    /// Object serialization
    /**
     * This method will save/load \p this into the archive \p ar.
     *
     * @param ar target archive.
     *
     * @throws unspecified any exception thrown by the serialization of the expression and of primitive types.
     */
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        // invoke serialization of the base class
        ar &boost::serialization::base_object<expression<T>>(*this);
        ar &m_weights;
        ar &m_weights_symbols;
    }

    // Delete ephemeral constants methods.
    void set_eph_val(const std::vector<T> &) = delete;
    void set_eph_symb(const std::vector<T> &) = delete;

private:
    // For numeric computations
    template <typename U, typename std::enable_if<std::is_same<U, double>::value || is_gdual<U>::value, int>::type = 0>
    U kernel_call(std::vector<U> &function_in, unsigned idx, unsigned node_id, unsigned weight_idx) const
    {
        // Weights (we transform the inputs a,b,c,d,e in w_1 a, w_2 b, w_3 c, etc...)
        for (auto j = 0u; j < this->_get_arity(node_id); ++j) {
            function_in[j] = function_in[j] * m_weights[weight_idx + j];
        }
        return this->get_f()[this->get()[idx]](function_in);
    }

    // For the symbolic expression
    template <typename U, typename std::enable_if<std::is_same<U, std::string>::value, int>::type = 0>
    U kernel_call(std::vector<U> &function_in, unsigned idx, unsigned node_id, unsigned weight_idx) const
    {
        // Weights (we transform the inputs x,y in (w_1*x), (w_2*y). The parenthesis is necessary
        // to avoid false representations such as w_1*x/w_2*y.
        for (auto j = 0u; j < this->_get_arity(node_id); ++j) {
            function_in[j] = "(" + m_weights_symbols[weight_idx + j] + "*" + function_in[j] + ")";
        }
        return this->get_f()[this->get()[idx]](function_in);
    }

    std::vector<T> m_weights;
    std::vector<std::string> m_weights_symbols;
};

} // end of namespace dcgp

#endif // DCGP_EXPRESSION_H
