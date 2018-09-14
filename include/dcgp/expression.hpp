#ifndef DCGP_EXPRESSION_H
#define DCGP_EXPRESSION_H

#include <algorithm>
#include <audi/audi.hpp>
#include <initializer_list>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <dcgp/kernel.hpp>
#include <dcgp/type_traits.hpp>

namespace dcgp
{

/// A dCGP expression
/**
 * This class represents a mathematical expression as encoded using CGP and
 * contains algorithms to compute its value (numerical and symbolical) and its
 * derivatives as well as to mutate the expression.
 *
 * @tparam T expression type. Can be double, or a gdual type.
 */
template <typename T>
class expression
{

private:
    // Static checks.
    static_assert(std::is_same<T, double>::value || is_gdual<T>::value,
                  "A d-CGP expression can only be operating on doubles or gduals");
    // SFINAE dust
    template <typename U>
    using functor_enabler = typename std::enable_if<
        std::is_same<U, double>::value || is_gdual<T>::value || std::is_same<U, std::string>::value, int>::type;

public:
    /// Constructor
    /** Constructs a dCGP expression with variable arity
     *
     * @param[in] n number of inputs (independent variables).
     * @param[in] m number of outputs (dependent variables).
     * @param[in] r number of rows of the dCGP.
     * @param[in] c number of columns of the dCGP.
     * @param[in] l number of levels-back allowed in the dCGP.
     * @param[in] arity arities of the basis functions for each column.
     * @param[in] f function set. An std::vector of dcgp::kernel<expression::type>.
     * @param[in] seed seed for the random number generator (initial expression
     * and mutations depend on this).
     */
    expression(unsigned n,                  // n. inputs
               unsigned m,                  // n. outputs
               unsigned r,                  // n. rows
               unsigned c,                  // n. columns
               unsigned l,                  // n. levels-back
               std::vector<unsigned> arity, // basis functions' arity
               std::vector<kernel<T>> f,    // functions
               unsigned seed                // seed for the pseudo-random numbers
               )
        : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_arity(arity), m_f(f), m_e(seed)
    {
        // Sanity checks
        sanity_checks();
        // Initializing bounds and chromosome
        init_bounds_and_chromosome();
        // We generate a random chromosome (expression)
        for (auto i = 0u; i < m_x.size(); ++i) {
            m_x[i] = std::uniform_int_distribution<unsigned>(m_lb[i], m_ub[i])(m_e);
        }
        update_data_structures();
    }

    /// Constructor
    /** Constructs a dCGP expression with uniform arity
     *
     * @param[in] n number of inputs (independent variables).
     * @param[in] m number of outputs (dependent variables).
     * @param[in] r number of rows of the dCGP.
     * @param[in] c number of columns of the dCGP.
     * @param[in] l number of levels-back allowed in the dCGP.
     * @param[in] arity arity of the basis functions.
     * @param[in] f function set. An std::vector of dcgp::kernel<expression::type>.
     * @param[in] seed seed for the random number generator (initial expression
     * and mutations depend on this).
     */
    expression(unsigned n,               // n. inputs
               unsigned m,               // n. outputs
               unsigned r,               // n. rows
               unsigned c,               // n. columns
               unsigned l,               // n. levels-back
               unsigned arity,           // basis functions' arity
               std::vector<kernel<T>> f, // functions
               unsigned seed             // seed for the pseudo-random numbers
               )
        : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_f(f), m_e(seed)
    {
        // We fill the arity vector with the same number (uniform arity)
        m_arity = std::vector<unsigned>(m_c, arity);
        // Sanity checks
        sanity_checks();
        // Initializing bounds and chromosome
        init_bounds_and_chromosome();
        // We generate a random chromosome (expression)
        for (auto i = 0u; i < m_x.size(); ++i) {
            m_x[i] = std::uniform_int_distribution<unsigned>(m_lb[i], m_ub[i])(m_e);
        }
        update_data_structures();
    }

    /// Virtual destructor (necessary for inheritance?)
    virtual ~expression(){};

    /// Evaluates the dCGP expression
    /**
     * This evaluates the dCGP expression. According to the template parameter
     * it will compute the value (double) the Taylor expansion (gdual) or a
     * symbolic representation (std::string). Any other type will result in a
     * compilation-time error (SFINAE).
     *
     * @param[point] an std::vector containing the values where the dCGP
     * expression has to be computed (doubles, gduals or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename U, functor_enabler<U> = 0>
    std::vector<U> operator()(const std::vector<U> &point) const
    {
        if (point.size() != m_n) {
            throw std::invalid_argument("Input size is incompatible");
        }
        std::vector<U> retval(m_m);
        std::vector<U> node(m_n + m_r * m_c);
        std::vector<U> function_in;

        for (auto i : m_active_nodes) {
            if (i < m_n) {
                node[i] = point[i];
            } else {
                unsigned arity = get_arity(i);
                function_in.resize(arity);
                unsigned idx = m_node_x_idx[i]; // position in the chromosome of the current node
                for (auto j = 0u; j < arity; ++j) {
                    function_in[j] = node[m_x[idx + j + 1u]];
                }
                node[i] = m_f[m_x[idx]](function_in);
            }
        }
        for (auto i = 0u; i < m_m; ++i) {
            retval[i] = node[m_x[m_x.size() - m_m + i]];
        }
        return retval;
    }

    /// Evaluates the dCGP expression
    /**
     * This evaluates the dCGP expression. According to the template parameter
     * it will compute the value (double) the Taylor expansion (gdual) or a
     * symbolic representation (std::string). Any other type will result in a
     * compilation-time error (SFINAE). This is identical to the other overload
     * and is provided only for convenience
     *
     * @param[in] in an initializer list containing the values where the dCGP
     * expression has to be computed (doubles, gduals or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename U, functor_enabler<U> = 0>
    std::vector<U> operator()(const std::initializer_list<U> &in) const
    {
        std::vector<U> dummy(in);
        return (*this)(dummy);
    }

    /// Sets the chromosome
    /** Sets a given chromosome as genotype for the expression and updates
     * the active nodes and active genes information accordingly
     *
     * @param[in] x the new cromosome
     *
     * @throw std::invalid_argument if the chromosome is incompatible with
     * the expression (n.inputs, n.outputs, levels-back, etc.)
     */
    void set(const std::vector<unsigned> &x)
    {
        if (!is_valid(x)) {
            throw std::invalid_argument("Chromosome is incompatible");
        }
        m_x = x;
        update_data_structures();
    }

    /// Sets the function gene of a node
    /** Sets for a valid node (i.e. not an input node) a new kernel
     *
     * @param[in] node_id the id of the node
     * @param[in] f_id the id of the kernel
     *
     * @throw std::invalid_argument if the *node_id* or *f_id* are invalid.
     */
    void set_f_gene(unsigned node_id, unsigned f_id)
    {
        if (f_id > m_f.size() - 1) {
            throw std::invalid_argument("You are trying to set a kernel id of: " + std::to_string(f_id)
                                        + ", but allowed values are [0 ... " + std::to_string(m_f.size() - 1)
                                        + "] since this CGP has " + std::to_string(m_f.size() - 1) + " kernels.");
        }
        if (node_id < m_n || node_id > m_n + m_c * m_r - 1u) {
            throw std::invalid_argument("You are trying to set the gene corresponding to a node_id: "
                                        + std::to_string(node_id) + ", but allowed values are [" + std::to_string(m_n)
                                        + " ... " + std::to_string(m_n + m_c * m_r - 1u) + "]");
        }
        auto gene_idx = m_node_x_idx[node_id];
        m_x[gene_idx] = f_id;
    }

    /// Gets the chromosome
    /**
     * Gets the chromosome encoding the current expression
     *
     * @return The chromosome
     */
    const std::vector<unsigned> &get() const
    {
        return m_x;
    }

    /// Gets the lower bounds
    /**
     * Gets the lower bounds for the genes
     *
     * @return An std::vector containing the lower bound for each gene
     */
    const std::vector<unsigned> &get_lb() const
    {
        return m_lb;
    }

    /// Gets the upper bounds
    /**
     * Gets the upper bounds for the genes
     *
     * @return An std::vector containing the upper bound for each gene
     */
    const std::vector<unsigned> &get_ub() const
    {
        return m_ub;
    }

    /// Gets the active genes
    /**
     * Gets the idx of the active genes in the current chromosome (numbering is
     * from 0)
     *
     * @return An std::vector containing the idx of the active genes in the
     * current chromosome
     */
    const std::vector<unsigned> &get_active_genes() const
    {
        return m_active_genes;
    }

    /// Gets the active nodes
    /**
     * Gets the idx of the active nodes in the current chromosome.
     * The numbering starts from 0 at the first input node to then follow PPSN
     * tutorial from Miller
     *
     * @return An std::vector containing the idx of the active nodes
     */
    const std::vector<unsigned> &get_active_nodes() const
    {
        return m_active_nodes;
    }

    /// Gets the number of inputs
    /**
     * Gets the number of inputs of the dCGP expression
     *
     * @return the number of inputs
     */
    unsigned get_n() const
    {
        return m_n;
    }

    /// Gets the number of outputs
    /**
     * Gets the number of outputs of the dCGP expression
     *
     * @return the number of outputs
     */
    unsigned get_m() const
    {
        return m_m;
    }

    /// Gets the number of rows
    /**
     * Gets the number of rows of the dCGP expression
     *
     * @return the number of rows
     */
    unsigned get_r() const
    {
        return m_r;
    }

    /// Gets the number of columns
    /**
     * Gets the number of columns of the dCGP expression
     *
     * @return the number of columns
     */
    unsigned get_c() const
    {
        return m_c;
    }

    /// Gets the number of levels-back
    /**
     * Gets the number of levels-back allowed for the dCGP expression
     *
     * @return the number of levels-back
     */
    unsigned get_l() const
    {
        return m_l;
    }

    /// Gets the arity
    /**
     * Gets the arity of the basis functions of the dCGP expression
     *
     * @return the arity
     */
    const std::vector<unsigned> &get_arity() const
    {
        return m_arity;
    }

    /// Gets the arity of a particular node
    /**
     * This does not perform any checks on the node_id. If an invalid node is queried ITS BAD an
     * invalid read will be performed.
     *
     * param[in] node_id id of the node
     * return the arity of that node
     *
     */
    unsigned get_arity(unsigned node_id) const
    {
        assert(node_id < m_r * m_c + m_n);
        unsigned col = (node_id - m_n) / m_r;
        return m_arity[col];
    }

    /// Gets the function set
    /**
     * Gets the set of functions used in the dCGP expression
     *
     * @return an std::vector of kernels
     */
    const std::vector<kernel<T>> &get_f() const
    {
        return m_f;
    }

    /// Gets node_x_idx
    /**
     * Gets node_x_idx, a vector containing the indexes in the chromosome where
     * nodes start expressing.
     *
     * @return an std::vector of kernels
     */
    const std::vector<unsigned> &get_node_x_idx() const
    {
        return m_node_x_idx;
    }

    /// Mutates randomly one gene
    /**
     * Mutates exactly one gene within its allowed bounds.
     *
     * @param[in] idx index of the gene to me mutated
     *
     * @throw std::invalid_argument if \p idx is too large
     */
    void mutate(unsigned idx)
    {
        if (idx >= m_x.size()) {
            throw std::invalid_argument("idx of gene to be mutated is out of bounds");
        }
        // If only one value is allowed for the gene, (lb==ub),
        // then we will not do anything as mutation does not apply
        if (m_lb[idx] < m_ub[idx]) {
            unsigned new_value;
            do {
                new_value = std::uniform_int_distribution<unsigned>(m_lb[idx], m_ub[idx])(m_e);
            } while (new_value == m_x[idx]);
            m_x[idx] = new_value;
            update_data_structures(); // TODO: unecessary if the gene is a function gene
        }
    }

    /// Mutates multiple genes at once
    /**
     * Mutates multiple genes within their allowed bounds.
     *
     * @param[in] idxs vector of indexes of the genes to me mutated
     *
     * @throw std::invalid_argument if \p idx is too large
     */
    void mutate(std::vector<unsigned> idxs)
    {
        bool flag = false;
        for (auto i = 0u; i < idxs.size(); ++i) {
            if (idxs[i] >= m_x.size()) {
                throw std::invalid_argument("idx of gene to be mutated is out of bounds");
            }
            // If only one value is allowed for the gene, (lb==ub),
            // then we will not do anything as mutation does not apply
            if (m_lb[idxs[i]] < m_ub[idxs[i]]) {
                unsigned new_value;
                do {
                    new_value = std::uniform_int_distribution<unsigned>(m_lb[idxs[i]], m_ub[idxs[i]])(m_e);
                } while (new_value == m_x[idxs[i]]);
                m_x[idxs[i]] = new_value;
                flag = true;
            }
        }
        if (flag) update_data_structures();
    }

    /// Mutates N random genes
    /**
     * Mutates a specified number of random genes within their bounds
     *
     * @param[in] N number of genes to be mutated
     *
     */
    void mutate_random(unsigned N)
    {
        bool flag = false;
        for (auto i = 0u; i < N; ++i) {
            // If only one value is allowed for the gene, (lb==ub),
            // then we will not do anything as mutation does not apply
            auto idx = std::uniform_int_distribution<unsigned>(0, m_lb.size() - 1)(m_e);
            if (m_lb[idx] < m_ub[idx]) {
                unsigned new_value;
                do {
                    new_value = std::uniform_int_distribution<unsigned>(m_lb[idx], m_ub[idx])(m_e);
                } while (new_value == m_x[idx]);
                m_x[idx] = new_value;
                flag = true;
            }
        }
        if (flag) update_data_structures();
    }

    /// Mutates active genes
    /**
     * Mutates \p N active genes within their allowed bounds.
     * The mutation can affect function genes, input genes and output genes.
     *
     * @param[in] N Number of active genes to be mutated
     *
     */
    void mutate_active(unsigned N = 1)
    {
        for (auto i = 0u; i < N; ++i) {
            unsigned idx
                = std::uniform_int_distribution<unsigned>(0, static_cast<unsigned>(m_active_genes.size() - 1u))(m_e);
            idx = m_active_genes[idx];
            mutate(idx);
        }
    }

    /// Mutates one of the active function genes
    /**
     * Mutates exactly one of the active function genes within its allowed bounds.
     */
    void mutate_active_fgene(unsigned N = 1u)
    {
        // If no active function gene exists, do nothing
        if (m_active_genes.size() > m_m) {
            for (auto i = 0u; i < N; ++i) {
                unsigned node_id = 0u;
                while (node_id < m_n) { // we get a random active node (there will be one that is not an input node)
                    node_id = m_active_nodes[std::uniform_int_distribution<unsigned>(
                        0, static_cast<unsigned>(m_active_nodes.size() - 1u))(m_e)];
                }
                // Since the first gene, for each node, is the function gene, we just mutate on that position
                mutate(m_node_x_idx[node_id]);
            }
        }
    }

    /// Mutates one of the active connection genes
    /**
     * Mutates exactly one of the active connection genes within its allowed
     * bounds.
     */
    void mutate_active_cgene(unsigned N = 1u)
    {
        // If no active function gene exists, do nothing
        if (m_active_genes.size() > m_m) {
            for (auto i = 0u; i < N; ++i) {
                unsigned idx = 0u;
                while (idx < m_n) { // we get a random active node (there will be one that is not an input node)
                    idx = m_active_nodes[std::uniform_int_distribution<unsigned>(
                        0, static_cast<unsigned>(m_active_nodes.size() - 1u))(m_e)];
                }
                idx = m_node_x_idx[idx]
                      + std::uniform_int_distribution<unsigned>(1, static_cast<unsigned>(get_arity(idx)))(m_e);
                mutate(idx);
            }
        }
    }

    /// Mutates one of the active output genes
    /**
     * Mutates exactly one of the output genes within its allowed bounds.
     */
    void mutate_ogene(unsigned N = 1)
    {
        unsigned idx;
        if (m_m > 1) {
            for (auto i = 0u; i < N; ++i) {
                idx = std::uniform_int_distribution<unsigned>(static_cast<unsigned>(m_active_genes.size() - m_m),
                                                              static_cast<unsigned>(m_active_genes.size() - 1u))(m_e);
            }

        } else {
            idx = static_cast<unsigned>(m_active_genes.size() - 1u);
        }
        idx = m_active_genes[idx];
        mutate(idx);
    }

    /// Checks if a given node is active
    /**
     *
     * @param[in] node_id the node to be checked
     *
     * @return True if the node *node_id* is active in the CGP expression.
     */
    bool is_active(const unsigned node_id) const
    {
        return (std::find(m_active_nodes.begin(), m_active_nodes.end(), node_id) != m_active_nodes.end());
    }

    /// Overloaded stream operator
    /**
     * Will return a formatted string containing a human readable representation
     * of the class
     *
     * @return std::string containing a human-readable representation of the
     * problem.
     */
    friend std::ostream &operator<<(std::ostream &os, const expression &d)
    {
        audi::stream(os, "d-CGP Expression:\n");
        audi::stream(os, "\tNumber of inputs:\t\t", d.m_n, '\n');
        audi::stream(os, "\tNumber of outputs:\t\t", d.m_m, '\n');
        audi::stream(os, "\tNumber of rows:\t\t\t", d.m_r, '\n');
        audi::stream(os, "\tNumber of columns:\t\t", d.m_c, '\n');
        audi::stream(os, "\tNumber of levels-back allowed:\t", d.m_l, '\n');
        audi::stream(os, "\tBasis function arity:\t\t", d.m_arity, '\n');
        audi::stream(os, "\tStart of the gene expressing the node:\t\t", d.m_node_x_idx, '\n');
        audi::stream(os, "\n\tResulting lower bounds:\t", d.m_lb);
        audi::stream(os, "\n\tResulting upper bounds:\t", d.m_ub, '\n');
        audi::stream(os, "\n\tCurrent expression (encoded):\t", d.m_x, '\n');
        audi::stream(os, "\tActive nodes:\t\t\t", d.m_active_nodes, '\n');
        audi::stream(os, "\tActive genes:\t\t\t", d.m_active_genes, '\n');
        audi::stream(os, "\n\tFunction set:\t\t\t", d.m_f, '\n');
        return os;
    }

    /// Validity of a chromosome
    /**
     * Checks if a chromosome (i.e. a sequence of integers) is a valid expression
     * by verifying its length and the bounds
     *
     * @param[in] x chromosome
     */
    bool is_valid(const std::vector<unsigned> &x) const
    {
        // Checking for length
        if (x.size() != m_lb.size()) {
            return false;
        }

        // Checking for bounds on all genes
        for (auto i = 0u; i < x.size(); ++i) {
            if ((x[i] > m_ub[i]) || (x[i] < m_lb[i])) {
                return false;
            }
        }
        return true;
    }

protected:
    /// Updates the class data that depend on the chromosome
    /**
     * Some of the expression data depend on the chromosome. This is the case, for example,
     * of the active nodes and active genes. Each time the chromosome is changed, these structures need also to be
     * changed. A call to this method takes care of this. In derived classes (such as for example expression_ann), one
     * can add more of these chromosome dependant data, and will thus need to override this method, making sure to still
     * have it called by the new method and adding there the new data book-keeping.
     */

    virtual void update_data_structures()
    {
        assert(m_x.size() == m_lb.size());

        // First we update the active nodes
        std::vector<unsigned> current(m_m), next;
        m_active_nodes.clear();

        // At the beginning, current contains only the node connected to the output nodes
        for (auto i = 0u; i < m_m; ++i) {
            current[i] = m_x[m_x.size() - m_m + i];
        }
        do {
            m_active_nodes.insert(m_active_nodes.end(), current.begin(), current.end());

            for (auto node_id : current) {
                if (node_id >= m_n) // we skip the input nodes as they do
                                    // not have any connection
                {
                    auto node_arity = get_arity(node_id);
                    for (auto i = 1u; i <= node_arity; ++i) {
                        next.push_back(m_x[m_node_x_idx[node_id] + i]);
                    }
                } else {
                    m_active_nodes.push_back(node_id);
                }
            }
            // We remove duplicates to avoid processing them and thus having a 2^N
            // complexity
            std::sort(next.begin(), next.end());
            next.erase(std::unique(next.begin(), next.end()), next.end());
            current = next;
            next.clear();
        } while (current.size() > 0);

        // We remove duplicates and keep m_active_nodes sorted
        std::sort(m_active_nodes.begin(), m_active_nodes.end());
        m_active_nodes.erase(std::unique(m_active_nodes.begin(), m_active_nodes.end()), m_active_nodes.end());

        // Then the active genes
        m_active_genes.clear();
        for (auto i = 0u; i < m_active_nodes.size(); ++i) {
            auto node_id = m_active_nodes[i];
            if (node_id >= m_n) {
                for (auto j = 0u; j <= get_arity(node_id); ++j) {
                    m_active_genes.push_back(m_node_x_idx[node_id] + j);
                }
            }
        }
        // Output genes are always active
        for (auto i = 0u; i < m_m; ++i) {
            m_active_genes.push_back(static_cast<unsigned>(m_x.size()) - m_m + i);
        }
    }

private:
    void sanity_checks()
    {
        if (m_n == 0) throw std::invalid_argument("Number of inputs is 0");
        if (m_m == 0) throw std::invalid_argument("Number of outputs is 0");
        if (m_c == 0) throw std::invalid_argument("Number of columns is 0");
        if (m_r == 0) throw std::invalid_argument("Number of rows is 0");
        if (m_l == 0) throw std::invalid_argument("Number of level-backs is 0");
        if (m_arity.size() != m_c)
            throw std::invalid_argument("The arity vector size (" + std::to_string(m_arity.size()) + ") must be the same as the number of columns (" + std::to_string(m_c) + ")");
        if (std::any_of(m_arity.begin(), m_arity.end(), [](unsigned a) { return a == 0; })) {
            throw std::invalid_argument("Basis functions arity cannot be zero");
        }
        if (m_f.size() == 0) throw std::invalid_argument("Number of basis functions is 0");
    }
    void init_bounds_and_chromosome()
    {
        // Chromosome size is r*c + sum(arity)*r + m
        auto size = m_r * m_c + m_r * std::accumulate(m_arity.begin(), m_arity.end(), 0u) + m_m;
        // Allocate bounds and chromosome and gene position
        m_x = std::vector<unsigned>(size, 0u);
        m_lb = std::vector<unsigned>(size, 0u);
        m_ub = std::vector<unsigned>(size, 0u);
        m_node_x_idx = std::vector<unsigned>(m_r * m_c + m_n, 0u);

        // We loop over all nodes and set function and connection genes
        unsigned k = 0u;
        for (auto i = 0u; i < m_c; ++i) {     // column first
            for (auto j = 0u; j < m_r; ++j) { // then rows
                // Function gene (lower bounds are all 0u)
                m_ub[k] = static_cast<unsigned>(m_f.size() - 1u);
                k++;
                // Connections genes
                for (auto l = 0u; l < m_arity[i]; ++l) {
                    m_ub[k] = m_n + i * m_r - 1u;
                    if (i >= m_l) { // only if level-backs allow a lower bound exists
                        m_lb[k] = m_n + m_r * (i - m_l);
                    }
                    k++;
                }
            }
        }
        // Bounds for the output genes
        for (auto i = size - m_m; i < size; ++i) {
            m_ub[i] = m_n + m_r * m_c - 1u;
            if (m_l <= m_c) {
                m_lb[i] = m_n + m_r * (m_c - m_l);
            }
        }
        // We compute the position of genes expressing a given node
        for (auto node_id = 0u; node_id < m_node_x_idx.size(); ++node_id) {
            if (node_id < m_n) {
                m_node_x_idx[node_id]
                    = 0u; // We put some unused values for the input nodes as they have no gene representation
            } else {
                unsigned col = (node_id - m_n) / m_r;
                unsigned row = (node_id - m_n) % m_r;
                unsigned acc = 0u;
                for (auto j = 0u; j < col; ++j) {
                    acc = acc + m_arity[j];
                }
                acc *= m_r;
                m_node_x_idx[node_id] = acc + row * m_arity[col] + (node_id - m_n);
            }
        }
    }

    // number of inputs
    unsigned m_n;
    // number of outputs
    unsigned m_m;
    // number of rows
    unsigned m_r;
    // number of columns
    unsigned m_c;
    // number of levels_back allowed
    unsigned m_l;
    // function arity
    std::vector<unsigned> m_arity;

    // the functions allowed
    std::vector<kernel<T>> m_f;
    // lower and upper bounds on all genes
    std::vector<unsigned> m_lb;
    std::vector<unsigned> m_ub;
    // active nodes idx (guaranteed to be always sorted)
    std::vector<unsigned> m_active_nodes;
    // active genes idx
    std::vector<unsigned> m_active_genes;
    // the encoded chromosome
    std::vector<unsigned> m_x;
    // The starting index in the chromosome of the genes expressing a node
    std::vector<unsigned> m_node_x_idx;
    // the random engine for the class
    std::default_random_engine m_e;
    // The expression type
    using type = T;
};

} // end of namespace dcgp

#endif // DCGP_EXPRESSION_H
