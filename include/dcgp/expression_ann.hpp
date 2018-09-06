#ifndef DCGP_EXPRESSION_ANN_H
#define DCGP_EXPRESSION_ANN_H

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

/// A dCGP-ANN expression
/**
 * This class represents an artificial neural network as a differentiable Cartesian Genetic
 * program. It add weights, biases and backward automated differentiation to the class
 * dcgp::expression.
 *
 * @tparam T expression type. Can only be a float type
 */
template <typename T>
class expression_ann : public expression<T>
{

private:
    template <typename U>
    using enable_double_string =
        typename std::enable_if<std::is_same<U, double>::value || std::is_same<U, std::string>::value, int>::type;
    template <typename U>
    using enable_double = typename std::enable_if<std::is_same<U, double>::value, int>::type;

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
    expression_ann(unsigned n,               // n. inputs
                   unsigned m,               // n. outputs
                   unsigned r,               // n. rows
                   unsigned c,               // n. columns
                   unsigned l,               // n. levels-back
                   unsigned arity,           // basis functions' arity
                   std::vector<kernel<T>> f, // functions
                   unsigned seed             // seed for the pseudo-random numbers
                   )
        : expression<T>(n, m, r, c, l, arity, f, seed), m_weights(r * c * arity, T(1.)), m_biases(r * c, T(0.))
    {
        for (const auto &ker : f) {
            if (ker.get_name() != "tanh" && ker.get_name() != "sig" && ker.get_name() != "ReLu") {
                throw std::invalid_argument("Only tanh, sig and ReLu Kernels are valid for dCGP-ANN expressions");
            }
        }
        for (auto i = 0u; i < r * c; ++i) {
            for (auto j = 0u; j < arity; ++j) {
                m_weights_symbols.push_back("w" + std::to_string(i + n) + "_" + std::to_string(j));
            }
        }
        for (auto i = 0u; i < r * c; ++i) {
            m_biases_symbols.push_back("b" + std::to_string(i + n));
        }
        update_data_structures();
    }

    /// Evaluates the dCGP-ANN expression
    /**
     * This evaluates the dCGP-ANN expression. According to the template parameter
     * it will compute the value (U) or a symbolic
     * representation (std::string). Any other type will result in a compilation-time
     * error (SFINAE).
     *
     * @param[point] in an std::vector containing the values where the dCGP-ANN expression has
     * to be computed (doubles or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename U, enable_double_string<U> = 0>
    std::vector<U> operator()(const std::vector<U> &point) const
    {
        std::vector<U> retval(this->get_m());
        auto node = fill_nodes(point);
        for (auto i = 0u; i < this->get_m(); ++i) {
            retval[i] = node[this->get()[(this->get_rows() * this->get_cols()) * (this->get_arity() + 1) + i]];
        }
        return retval;
    }

    /// Evaluates the  dCGP-ANN  expression
    /**
     * This evaluates the dCGP-ANN expression. According to the template parameter
     * it will compute the value (U) or a symbolic
     * representation (std::string). Any other type will result in a compilation-time
     * error (SFINAE).
     *
     * @param[point] in is an initializer list containing the values where the dCGP expression has
     * to be computed (U, or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename U, enable_double_string<U> = 0>
    std::vector<U> operator()(const std::initializer_list<U> &point) const
    {
        std::vector<U> dummy(point);
        return (*this)(dummy);
    }

    /// Evaluates the mean square error and its gradient
    /**
     * Returns the mean squared error and its gradient with respect to weights and biases.
     *
     * @param[point] The input data (single point)
     * @param[prediction] The predicted output (single point)
     * @return the mse, the gradient of the mse w.r.t. all weights (also inactive) and the gradient of the mse w.r.t all
     * biases
     */
    template <typename U, enable_double<U> = 0>
    std::tuple<U, std::vector<U>, std::vector<U>> mse(const std::vector<U> &point, const std::vector<U> &prediction)
    {
        if (point.size() != this->get_n()) {
            throw std::invalid_argument("When computing the mse the point dimension (input) seemed wrong, it was: "
                                        + std::to_string(point.size())
                                        + " while I expected: " + std::to_string(this->get_n()));
        }
        if (prediction.size() != this->get_m()) {
            throw std::invalid_argument(
                "When computing the mse the prediction dimension (output) seemed wrong, it was: "
                + std::to_string(prediction.size()) + " while I expected: " + std::to_string(this->get_m()));
        }
        U value(U(0.));
        std::vector<U> gweights(m_weights.size(), U(0.));
        std::vector<U> gbiases(m_biases.size(), U(0.));

        // ------------------------------------------ Forward pass (takes roughly half of the time) --------------------
        // All active nodes outputs get computed as well as
        // the activation function derivatives
        auto n_nodes = this->get_n() + this->get_rows() * this->get_cols();
        std::vector<U> node(n_nodes, 0.), d_node(n_nodes, 0.);
        fill_nodes(point, node, d_node); // here is where the computations happen.

        // We add to node and node_d some virtual nodes containing the output values (x_i-\hat x_i) ^ 2, and its
        // derivative 2(x_i - \hat x_i) so that the acyclic graph is now computing all blocks of the mse (a dCGP-ANN
        // only computes outputs)
        for (decltype(this->get_m()) i = 0u; i < this->get_m(); ++i) {
            auto node_idx = this->get()[this->get().size() - this->get_m() + i];
            auto dummy = (node[node_idx] - prediction[i]);
            node.push_back(dummy * dummy);
            d_node.push_back(2 * dummy);
            value += dummy * dummy;
        }

        // ------------------------------------------ Backward pass (takes roughly the remaining half) -----------------
        // We iterate backward on all the active nodes (except the input nodes)
        // filling up the gradient information at each node for the incoming weights and relative bias
        for (auto it = this->get_active_nodes().rbegin(); it != this->get_active_nodes().rend(); ++it) {
            if (*it < this->get_n()) continue;
            // index of the node in the bias vector
            auto b_idx = *it - this->get_n();
            // index of the node in the weight vector
            auto w_idx = this->get_arity() * (*it - this->get_n());
            // index of the node in the chromosome
            auto c_idx = (*it - this->get_n()) * (this->get_arity() + 1);
            // index in the node/d_node vectors
            auto n_idx = *it;
            // We update the d_node information
            U cum = 0.;
            for (auto i = 0u; i < m_connected[*it].size(); ++i) {
                // If the node is not "virtual", that is not one of the m virtual nodes we added computing (x-x_i)^2
                if (m_connected[*it][i].first < this->get_n() + this->get_rows() * this->get_cols()) {
                    cum += m_weights[m_connected[*it][i].second] * d_node[m_connected[*it][i].first];
                } else {
                    auto n_out = m_connected[*it][i].first - (this->get_n() + this->get_rows() * this->get_cols());
                    cum += d_node[d_node.size() - this->get_m() + n_out];
                }
            }
            d_node[n_idx] *= cum;

            // fill gradients for weights and biases info
            for (auto i = 0u; i < this->get_arity(); ++i) {
                gweights[w_idx + i] = d_node[n_idx] * node[this->get()[c_idx + 1 + i]];
            }
            gbiases[b_idx] = d_node[n_idx];
        }

        return std::make_tuple(std::move(value), std::move(gweights), std::move(gbiases));
    }

    /// Evaluates the mean square error and its gradient
    /**
     * Returns the mean squared error and its gradient with respect to weights and biases.
     *
     * @param[points] The input data (a batch).
     * @param[predictions] The predicted outputs (a batch).
     * @return the mse, the gradient of the mse w.r.t. all weights (also inactive) and the gradient of the mse w.r.t all
     * biases.
     */
    template <typename U, enable_double<U> = 0>
    std::tuple<U, std::vector<U>, std::vector<U>> mse(const std::vector<std::vector<U>> &points,
                                                      const std::vector<std::vector<U>> &predictions)
    {
        if (points.size() != predictions.size()) {
            throw std::invalid_argument("Data and label size mismatch data size is: " + std::to_string(points.size())
                                        + " while label size is: " + std::to_string(predictions.size()));
        }
        if (points.size() == 0) {
            throw std::invalid_argument("Data size cannot be zero");
        }
        return mse<U>(points.begin(), points.end(), predictions.begin());
    }

    /// Stochastic gradient descent
    /**
     * Performs one "epoch" of stochastic gradient descent using mean square error
     *
     * @param[points] The input data (a batch).
     * @param[predictions] The predicted outputs (a batch).
     * @param[l_rate] The learning rate.
     * @param[batch_size] The batch size.
     *
     * @throws std::invalid_argument if the *data* and *label* size do not match or is zero, or if *l_rate* is not
     * positive.
     */
    template <typename U, enable_double<U> = 0>
    void sgd(const std::vector<std::vector<U>> &points, const std::vector<std::vector<U>> &predictions, double l_rate,
             unsigned batch_size)
    {
        if (points.size() != predictions.size()) {
            throw std::invalid_argument("Data and label size mismatch data size is: " + std::to_string(points.size())
                                        + " while label size is: " + std::to_string(predictions.size()));
        }
        if (points.size() == 0) {
            throw std::invalid_argument("Data size cannot be zero");
        }
        if (l_rate <= 0) {
            throw std::invalid_argument("The learning rate must be a positive number, while: " + std::to_string(l_rate)
                                        + " was detected.");
        }
        auto dfirst = points.begin();
        auto dlast = points.end();
        auto lfirst = predictions.begin();
        while (dfirst != dlast) {
            if (dfirst + batch_size > dlast) {
                update_weights<U>(dfirst, dlast, lfirst, l_rate);
                dfirst = dlast;
            } else {
                update_weights<U>(dfirst, dfirst + batch_size, lfirst, l_rate);
                dfirst += batch_size;
                lfirst += batch_size;
            }
        }
    }

    /// Overloaded stream operator
    /**
     * Will return a formatted string containing a human readable representation
     * of the class
     *
     * @return std::string containing a human-readable representation of the problem.
     */
    friend std::ostream &operator<<(std::ostream &os, const expression_ann &d)
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
        audi::stream(os, "\tBiases:\t\t\t\t", d.m_biases, '\n');

        return os;
    }

    /**
     * \defgroup Managing Weights and Biases
     */
    /*@{*/

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

    /// Sets a weight
    /**
     * Sets a connection weight to a new value
     *
     * @param[idx] index of the weight to be changed.
     * @param[w] value of the weight to be changed.
     *
     * @throws std::invalid_argument if the node_id or input_id are not valid
     */
    void set_weight(typename std::vector<T>::size_type idx, const T &w)
    {
        m_weights[idx] = w;
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
    T get_weight(typename std::vector<T>::size_type node_id, typename std::vector<T>::size_type input_id) const
    {
        if (node_id < this->get_n() || node_id >= this->get_n() + this->get_rows() * this->get_cols()) {
            throw std::invalid_argument(
                "Requested node id does not exist or does not have a weight (e.g. input nodes)");
        }
        if (input_id >= this->get_arity()) {
            throw std::invalid_argument("Requested input exceeds the function arity");
        }
        auto idx = (node_id - this->get_n()) * this->get_arity() + input_id;
        return m_weights[idx];
    }

    /// Gets a weight
    /**
     * Gets the value of a connection weight
     *
     * @param[in] idx index of the weight
     *
     */
    T get_weight(typename std::vector<T>::size_type idx) const
    {
        return m_weights[idx];
    }

    /// Gets the weights
    /**
     * Gets the values of all the weights.
     *
     * @return an std::vector containing all the weights
     */
    const std::vector<T> &get_weights() const
    {
        return m_weights;
    }

    /// Randomises all weights
    /**
     * Set all weights to a normally distributed number
     *
     * @param[mean] the mean of the normal distribution.
     * @param[std] the standard deviation of the normal distribution.
     * @param[seed] the seed to generate the new weights (by default its randomly generated).
     *
     */
#if !defined(DCGP_DOXYGEN_INVOKED)
    void randomise_weights(double mean = 0, double std = 0.1,
                           std::random_device::result_type seed = std::random_device{}())
    {
        std::mt19937 gen{seed};
        std::normal_distribution<T> nd{mean, std};
        for (auto &w : m_weights) {
            w = nd(gen);
        }
    }
#else
    void randomise_weights(double mean = 0, double std = 0.1, std::random_device::result_type seed = random_number) {}
#endif

    /// Sets a bias
    /**
     * Sets a node bias to a new value
     *
     * @param[idx] index of the bias to be changed.
     * @param[w] value of the new bias.
     *
     */
    void set_bias(typename std::vector<T>::size_type idx, const T &w)
    {
        m_biases[idx] = w;
    }

    /// Sets all biases
    /**
     * Sets all the nodes biases at once
     *
     * @param[in] bs an std::vector containing all the biases to set
     *
     * @throws std::invalid_argument if the input vector dimension is not valid (r*c)
     */
    void set_biases(const std::vector<T> &bs)
    {
        if (bs.size() != m_biases.size()) {
            throw std::invalid_argument("The vector of biases has the wrong dimension");
        }
        m_biases = bs;
    }

    /// Gets a bias
    /**
     * Gets the value of a bias
     *
     * @param[in] idx index of the bias
     *
     */
    T get_bias(typename std::vector<T>::size_type idx) const
    {
        return m_biases[idx];
    }

    /// Gets the biases
    /**
     * Gets the values of all the biases
     *
     * @return an std::vector containing all the biases
     */
    const std::vector<T> &get_biases() const
    {
        return m_biases;
    }

    /// Randomises all biases
    /**
     * Set all biases to a normally distributed number
     *
     * @param[in] mean the mean of the normal distribution
     * @param[in] std the standard deviation of the normal distribution
     * @param[in] seed the seed to generate the new biases (by default its randomly generated)
     *
     */
#if !defined(DCGP_DOXYGEN_INVOKED)
    void randomise_biases(double mean = 0., double std = 0.1,
                          std::random_device::result_type seed = std::random_device{}())
    {
        std::mt19937 gen{seed};
        std::normal_distribution<T> nd{mean, std};
        for (auto &b : m_biases) {
            b = nd(gen);
        }
    }
#else
    void randomise_biases(double mean = 0, double std = 0.1, std::random_device::result_type seed = random_number) {}
#endif

    /*@}*/

private:
    // For numeric computations
    template <typename U, enable_double<U> = 0>
    U kernel_call(std::vector<U> &function_in, unsigned idx, unsigned weight_idx, unsigned bias_idx) const
    {
        // Weights (we transform the inputs a,b,c,d,e in w_1 a, w_2 b, w_3 c, etc...)
        for (auto j = 0u; j < this->get_arity(); ++j) {
            function_in[j] = function_in[j] * m_weights[weight_idx + j];
        }
        // Biases (we add to the first input a bias so that a,b,c,d,e goes in c, etc...))
        function_in[0] += m_biases[bias_idx];
        // We compute the node function that will, for example, map w_1 a + bias, w_2 b, w_3 c,... into f(w_1 a +
        // w_2 b + w_3 c + ... + bias)
        return this->get_f()[this->get()[idx]](function_in);
    }

    // For the symbolic expression
    template <typename U, typename std::enable_if<std::is_same<U, std::string>::value, int>::type = 0>
    U kernel_call(std::vector<U> &function_in, unsigned idx, unsigned weight_idx, unsigned bias_idx) const
    {
        // Weights
        for (auto j = 0u; j < this->get_arity(); ++j) {
            function_in[j] = m_weights_symbols[weight_idx + j] + "*" + function_in[j];
        }
        // Biases
        function_in[0] = m_biases_symbols[bias_idx] + "+" + function_in[0];
        return this->get_f()[this->get()[idx]](function_in);
    }

    // computes node to evaluate the expression
    template <typename U, enable_double_string<U> = 0>
    std::vector<U> fill_nodes(const std::vector<U> &in) const
    {
        if (in.size() != this->get_n()) {
            throw std::invalid_argument("Input size is incompatible");
        }
        std::vector<U> node(this->get_n() + this->get_rows() * this->get_cols());
        std::vector<U> function_in(this->get_arity());
        for (auto i : this->get_active_nodes()) {
            if (i < this->get_n()) {
                node[i] = in[i];
            } else {
                // position in the chromosome of the current node
                unsigned idx = (i - this->get_n()) * (this->get_arity() + 1);
                // starting position in m_weights of the weights relative to the node
                unsigned weight_idx = (i - this->get_n()) * this->get_arity();
                // starting position in m_biases of the node bias
                unsigned bias_idx = i - this->get_n();
                for (auto j = 0u; j < this->get_arity(); ++j) {
                    function_in[j] = node[this->get()[idx + j + 1]];
                }
                node[i] = kernel_call(function_in, idx, weight_idx, bias_idx);
            }
        }
        return node;
    }

    // computes node and node_d to start backprop
    template <typename U, enable_double<U> = 0>
    void fill_nodes(const std::vector<U> &in, std::vector<U> &node, std::vector<U> &d_node) const
    {
        if (in.size() != this->get_n()) {
            throw std::invalid_argument("Input size is incompatible");
        }
        // Start
        std::vector<U> function_in(this->get_arity());
        for (auto i : this->get_active_nodes()) {
            if (i < this->get_n()) {
                node[i] = in[i];
                // We need d_node to have the same structure of node, hence we also
                // put some bogus entries fot the input nodes that actually do not have an activation function
                // hence no need/use/meaning for a derivative
                d_node[i] = 0.;
            } else {
                // position in the chromosome of the current node
                unsigned idx = (i - this->get_n()) * (this->get_arity() + 1);
                // starting position in m_weights of the weights relative to the node
                unsigned weight_idx = (i - this->get_n()) * this->get_arity();
                // starting position in m_biases of the node bias
                unsigned bias_idx = i - this->get_n();
                for (auto j = 0u; j < this->get_arity(); ++j) {
                    function_in[j] = node[this->get()[idx + j + 1]];
                }
                node[i] = kernel_call(function_in, idx, weight_idx, bias_idx);
                // take cares of d_node
                // sigmoid derivative is sig(1-sig)
                if (this->get_f()[this->get()[idx]].get_name() == "sig") {
                    d_node[i] = node[i] * (1. - node[i]);
                    // tanh derivative is 1 - tanh**2
                } else if (this->get_f()[this->get()[idx]].get_name() == "tanh") {
                    d_node[i] = 1. - node[i] * node[i];
                    // Relu derivative is 0 if relu<0, 1 otherwise
                } else if (this->get_f()[this->get()[idx]].get_name() == "ReLu") {
                    d_node[i] = (node[i] > 0.) ? 1. : 0.;
                }
            }
        }
        return;
    }

    // This overrides the base class update_data_structures and updates also the m_connected (as well as
    // m_active_nodes and genes). It is called upon construction and each time active genes are changed.
    void update_data_structures()
    {
        expression<T>::update_data_structures();
        m_connected.clear();
        m_connected.resize(this->get_n() + this->get_m() + this->get_rows() * this->get_cols());
        for (auto node_id : this->get_active_nodes()) {
            if (node_id >= this->get_n()) { // not for input nodes
                // position in the chromosome of the current node
                unsigned idx = (node_id - this->get_n()) * (this->get_arity() + 1);
                // loop over the genes representing connections
                for (auto i = idx + 1; i < idx + 1 + this->get_arity(); ++i) {
                    if (this->is_active(this->get()[i])) {
                        m_connected[this->get()[i]].push_back(
                            {node_id, (node_id - this->get_n()) * this->get_arity() + i - idx - 1});
                    }
                }
            }
        }
        // We now add the output nodes with ids starting from n + r * c. In this case the weight is not
        // relevant, hence we use the arbitrary value 0u as index in the weight vector.
        for (auto i = 0u; i < this->get_m(); ++i) {
            auto virtual_idx = this->get_n() + this->get_rows() * this->get_cols() + i;
            auto node_idx = this->get()[this->get().size() - this->get_m() + i];
            m_connected[node_idx].push_back({virtual_idx, 0u});
        }
    }

    /// Performs one weight/bias update
    /**
     * Updates m_weights and m_biases using gradient descent
     *
     * @param[dfirst] Start range for the data
     * @param[dlast] End range for the data
     * @param[lfirst] Start range for the labels
     * @param[lr] The learning rate
     *
     * @throws std::invalid_argument if the *data* and *label* size do not match or are zero, or if *lr* is not
     * positive.
     */
    template <typename U, enable_double<U> = 0>
    void update_weights(typename std::vector<std::vector<U>>::const_iterator dfirst,
                        typename std::vector<std::vector<U>>::const_iterator dlast,
                        typename std::vector<std::vector<U>>::const_iterator lfirst, U lr)
    {
        U coeff(lr / static_cast<U>(dlast - dfirst));
        while (dfirst != dlast) {
            auto mse_out = mse(*dfirst++, *lfirst++);
            std::transform(m_weights.begin(), m_weights.end(), std::get<1>(mse_out).begin(), m_weights.begin(),
                           [coeff](U a, U b) { return a - coeff * b; });
            std::transform(m_biases.begin(), m_biases.end(), std::get<2>(mse_out).begin(), m_biases.begin(),
                           [coeff](U a, U b) { return a - coeff * b; });
        }
    }

    template <typename U, enable_double<U> = 0>
    std::tuple<U, std::vector<U>, std::vector<U>> mse(typename std::vector<std::vector<U>>::const_iterator dfirst,
                                                      typename std::vector<std::vector<U>>::const_iterator dlast,
                                                      typename std::vector<std::vector<U>>::const_iterator lfirst)
    {
        U value(U(0.));
        std::vector<U> gweights(m_weights.size(), U(0.));
        std::vector<U> gbiases(m_biases.size(), U(0.));
        U dim = static_cast<U>(dlast - dfirst);

        while (dfirst != dlast) {
            auto mse_out = mse(*dfirst++, *lfirst++);
            value += std::get<0>(mse_out);
            std::transform(gweights.begin(), gweights.end(), std::get<1>(mse_out).begin(), gweights.begin(),
                           [dim](U a, U b) { return a + b / dim; });
            std::transform(gbiases.begin(), gbiases.end(), std::get<2>(mse_out).begin(), gbiases.begin(),
                           [dim](U a, U b) { return a + b / dim; });
        }
        value /= static_cast<U>(dim);

        return std::make_tuple(std::move(value), std::move(gweights), std::move(gbiases));
    }

private:
    std::vector<T> m_weights;
    std::vector<std::string> m_weights_symbols;

    std::vector<T> m_biases;
    std::vector<std::string> m_biases_symbols;

    // In order to be able to perform backpropagation on the dCGPANN program, we need to add
    // to the usual CGP data structures one that contains for each node the list of nodes
    // (and weights) it feeds into. We also need to add some virtual nodes (output nodes) 
    // computing the error components (x_i-\hat x_i) as to be able to get the mse deirvatives
    // Assigned virtual ids starting from n + r * c
    std::vector<std::vector<std::pair<unsigned, unsigned>>> m_connected;
}; // namespace dcgp

} // end of namespace dcgp

#endif // DCGP_EXPRESSION_H
