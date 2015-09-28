#ifndef DCGP_EXPRESSION_H
#define DCGP_EXPRESSION_H

#include <vector>
#include <string>
#include <map>
#include <random>
#include <initializer_list>
#include <stdexcept>
#include <audi/gdual.hpp>
#include <iostream>
#include <sstream>

#include "basis_function.hpp"
#include "std_overloads.hpp"


namespace dcgp {

/// A d-CGP expression
/**
 * This class represent a mathematical expression as encoded using CGP and contains
 * algorithms that compute its value (numerical and symbolical) and its derivatives 
 * as well as mutate the expression. 
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class expression {

private:
    // SFINAE dust
    template <typename T>
    using functor_enabler = typename std::enable_if<std::is_same<T,double>::value || std::is_same<T,gdual>::value || std::is_same<T,std::string>::value,int>::type;
public:   
    /// Constructor
    /** Constructs a d-CGP expression
     *
     * \param[in] n number of inputs (independent variables)
     * \param[in] m number of outputs (dependent variables)
     * \param[in] r number of rows of the cartesian cgp
     * \param[in] c number of columns of the cartesian cgp
     * \param[in] l number of levels-back allowed for the cartesian cgp
     * \param[in] f function set. An std::vector of dcgp::basis_function
     * \param[in] seed seed for the random number generator (initial expression  and mutations depend on this)
     */          
    expression(unsigned int n,                  // n. inputs
               unsigned int m,                  // n. outputs
               unsigned int r,                  // n. rows
               unsigned int c,                  // n. columns
               unsigned int l,                  // n. levels-back
               std::vector<basis_function> f,   // functions
               unsigned int seed                // seed for the pseudo-random numbers
               ) : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_f(f), m_lb((3 * m_r * m_c) + m_m, 0), m_ub((3 * m_r * m_c) + m_m, 0), m_x((3 * m_r * m_c) + m_m, 0), m_e(seed)
    {
        // Sanity checks
        if (n == 0) throw std::invalid_argument("Number of inputs is 0");
        if (m == 0) throw std::invalid_argument("Number of outputs is 0");
        if (c == 0) throw std::invalid_argument("Number of columns is 0");
        if (r == 0) throw std::invalid_argument("Number of rows is 0");
        if (l == 0) throw std::invalid_argument("Number of level-backs is 0");
        if (f.size()==0) throw std::invalid_argument("Number of basis functions is 0");

        // Bounds for the function genes
        for (auto i = 0u; i < (3 * m_r * m_c); i+=3) {
            m_ub[i] = f.size() - 1;
        }

        // Bounds for the output genes
        for (auto i = 3u * m_r * m_c; i < m_ub.size(); ++i) {
            m_ub[i] = m_n + m_r * m_c - 1;
            if (m_l <= m_c) {
                m_lb[i] = m_n + m_r * (m_c - m_l);
            }
        }

        // Bounds for the node connection genes 
        for (auto i = 0u; i < m_c; ++i) {
            for (auto j = 0u; j < m_r; ++j) {
                m_ub[((i * m_r) + j) * 3 + 1] = m_n + i * m_r - 1;
                m_ub[((i * m_r) + j) * 3 + 2] = m_n + i * m_r - 1;
                if (i >= m_l) {
                    m_lb[((i * m_r) + j) * 3 + 1] = m_n + m_r * (i - m_l);
                    m_lb[((i * m_r) + j) * 3 + 2] = m_n + m_r * (i - m_l);
                }
            }
        }

        // We generate a random expression
        for (auto i = 0u; i < m_x.size(); ++i)
        {
            m_x[i] = std::uniform_int_distribution<unsigned int>(m_lb[i], m_ub[i])(m_e);
        }
        update_active();
    }

    /// Sets the chromosome
    /** Sets a new chromosome as genotype for the expression and updates
     * the active nodes and active genes information
     *
     * \param[in] x The new cromosome
     *
     * @throw std::invalid_argument if the chromosome is incompatible with
     * the expression (n.inputs, n.outputs, levels-back, etc.)
     */
    void set(const std::vector<unsigned int> &x)
    {
        if(!is_valid(x))
        {
            throw std::invalid_argument("Chromosome is incompatible");
        }
        m_x = x;
        update_active();
    }

    /// Gets the chromosome
    /** 
     * Gets the chromosome encoding the current expression
     *
     * @return The chromosome
    */
    const std::vector<unsigned int> & get() const {return m_x;}

    /// Gets the active genes
    /** 
     * Gets the idx of the active genes in the current chromosome (numbering is from 0)
     *
     * @return An std::vector containing the idx of the active genes in the current chromosome
    */
    const std::vector<unsigned int> & get_active_genes() const {return m_active_genes;}
 
    /// Gets the active nodes
    /** 
     * Gets the idx of the active nodes in the current chromosome.
     * The numbering starts from 0 at the first input node to then follow PPSN tutorial from Miller
     *
     * @return An std::vector containing the idx of the active nodes
    */
    const std::vector<unsigned int> & get_active_nodes() const {return m_active_nodes;}
    
    /// Gets the number of inputs
    /** 
     * Gets the number of inputs of the c_CGP expression
     *
     * @return the number of inputs
    */
    unsigned int get_n() const {return m_n;}

    /// Gets the number of outputs
    /** 
     * Gets the number of outputs of the c_CGP expression
     *
     * @return the number of outputs
    */
    unsigned int get_m() const {return m_m;}

    /// Gets the functions
    /** 
     * Gets the functions
     *
     * @return an std::vector<basis_function>
    */
    const std::vector<basis_function>& get_f() const {return m_f;}

    /// Mutates one of the active genes
    /** 
     * Mutates exactly one of the active genes within its allowed bounds. The mutation
     * can affect a function gene or an input gene or an output gene.
     */
    void mutate_active()
    {
        unsigned int idx = std::uniform_int_distribution<unsigned int>(0, m_active_genes.size() - 1)(m_e);
        idx = m_active_genes[idx];

        if (m_lb[idx]<m_ub[idx]) // if only one value is allowed for the gene, then we will not do anything as mutation does not apply
        {
            unsigned int new_value = UINT_MAX;
            do 
            {
                new_value = std::uniform_int_distribution<unsigned int>(m_lb[idx], m_ub[idx])(m_e);
            } while (new_value == m_x[idx]);
            m_x[idx] = new_value;
            update_active();
        }
    }

    /// Evaluates the d-CGP expression
    /*
     * This evaluates the d-CGP expression. According to the template parameter
     * it will compute the value (double) the Taylor expansion (gdual) or a symbolic
     * representation (std::string). Any other type will result in a compilation-time
     * error (SFINAE).
     *
     * @param[in] in an std::vector containing the values where the d-CGP expression has
     * to be computed (doubles, gduals or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename T, functor_enabler<T> = 0>
    std::vector<T> operator()(const std::vector<T>& in) const
    {  
        if(in.size() != m_n)
        {
            throw std::invalid_argument("Input size is incompatible");
        }
        std::vector<T> retval(m_m);
        std::map<unsigned int, T> node;
        for (auto i : m_active_nodes) {
            if (i < m_n) 
            {
                node[i] = in[i];
            } else {
                unsigned int idx = (i - m_n) * 3;
                node[i] = m_f[m_x[idx]](node[m_x[idx + 1]], node[m_x[idx + 2]]);
            }
        }
        for (auto i = 0u; i<m_m; ++i)
        {
            retval[i] = node[m_x[(m_r * m_c) * 3 + i]];
        }
        return retval;
    }

    /// Evaluates the d-CGP expression
    /*
     * This evaluates the d-CGP expression. According to the template parameter
     * it will compute the value (double) the Taylor expansion (gdual) or a symbolic
     * representation (std::string). Any other type will result in a compilation-time
     * error (SFINAE). This is identical to the other overload and is provided only
     * for convenience
     *
     * @param[in] in an initializer list containing the values where the d-CGP expression has
     * to be computed (doubles, gduals or strings)
     *
     * @return The value of the function (an std::vector)
     */
    template <typename T, functor_enabler<T> = 0>
    std::vector<T> operator()(const std::initializer_list<T>& in) const
    {  
        std::vector<T> dummy(in);
        return (*this)(dummy);
    }


    /// Computes all derivatives of the expression up to order
    /** 
     * Using audi::gdual, this method returns the Taylor expansion of the 
     * d-CGP expression around the point \p in up to order \p order
     *
     * \param[in] in expansion point
     * \param[in] order maximum order for the Taylor expansion (and hence the derivatives)
     *
     * @return an std::vector<audi::gdual> containing the Taylor expansion 
     * derivatives can be extracted with retval.get_derivative({i,j,k, ... })
     *
     * @throw std::invalid_argument if the size of /p in is inconsistent with
     * the d-CGP number of inputs.
     */
    std::vector<audi::gdual> derivatives(const std::vector<double>& in, unsigned int order) const
    {  
        // We perform sanity checks
        if(in.size() != m_n)
        {
            throw std::invalid_argument("Input size is incompatible");
        }
        // We define the initial variables as xi = xi + dxi (dxi symbolic)
        std::vector<audi::gdual> in_expansion;
        for (auto i = 0u; i < in.size(); ++i) {
            in_expansion.emplace_back(in[i], "x" + std::to_string(i), (int)order);
        } 

        // We compute the CGP expression using gduals
        std::vector<audi::gdual> retval = (*this)(in_expansion);

        // If needed, we force the symbol set to contain all dxi
        // Needed as the CGP expression could be not using one of the inputs
        std::vector<std::string> symbols;
        for  (auto &expression : retval) {
            if (expression.get_symbol_set_size() != m_n) {
                if (symbols.size() == 0) {
                    for (auto i = 0u; i < in.size(); ++i) {
                        symbols.emplace_back("dx"+std::to_string(i));
                    } 
                }
                expression.extend_symbol_set(symbols);
            }
        }

        // We return the result
        return retval;
    }

    /// Overloaded stream operator
    /**
     * Will return a formatted string containing a human readable representation
     * of the class
     *
     * @return std::string containing a human-readable representation of the problem.
     */
    friend std::ostream &operator<<(std::ostream &s, const expression &d)
    {
        s << "d-CGP Expression:\n";
        s << "\tNumber of inputs:\t\t" << d.m_n << '\n';
        s << "\tNumber of outputs:\t\t" << d.m_m << '\n';
        s << "\tNumber of rows:\t\t\t" << d.m_r << '\n';
        s << "\tNumber of columns:\t\t" << d.m_c << '\n';
        s << "\tNumber of levels-back allowed:\t" << d.m_l << '\n';
        s << "\n\tResulting lower bounds:\t" << d.m_lb;
        s << "\n\tResulting upper bounds:\t" << d.m_ub << '\n';
        s << "\n\tCurrent expression (encoded):\t" << d.m_x << '\n';
        s << "\tActive nodes:\t\t\t" << d.m_active_nodes << '\n';
        s << "\tActive genes:\t\t\t" << d.m_active_genes << '\n';
        s << "\n\tFunction set:\t\t\t" << d.m_f << '\n';
        return s;
    }

protected:
    /// Validity of a chromosome
    /** 
     * Checks if a chromosome (i.e. a sequence of integers) is a valid expression
     * by checking its length and the bounds
     *
     * \param[in] x chromosome 
     */
    bool is_valid(const std::vector<unsigned int>& x) const
    {
        // Checking for length
        if (x.size() != m_lb.size()) {
            return false;
        }

        // Checking for bounds on all cenes
        for (auto i = 0u; i < x.size(); ++i) {
            if ((x[i] > m_ub[i]) || (x[i] < m_lb[i])) {
                return false;
            }
        }
        return true;
    }

    void update_active()
    {
        assert(m_x.size() == m_lb.size());

        // First we update the active nodes
        std::vector<unsigned int> current(m_m), next;
        m_active_nodes.clear();
        // At the beginning current contains only the output nodes connections
        for (auto i = 0u; i < m_m; ++i) {
            current[i] = m_x[3 * m_r * m_c + i];
        }
        do
        {
            m_active_nodes.insert(m_active_nodes.end(), current.begin(), current.end());

            for (auto node_id : current)
            {
                if (node_id >=m_n) // we insert the input nodes connections as they do not have any
                {
                    next.push_back(m_x[(node_id - m_n) * 3 + 1]);
                    next.push_back(m_x[(node_id - m_n) * 3 + 2]);
                }
                else{
                    m_active_nodes.push_back(node_id);
                }
            }
            // We remove duplicates to avoid processng them and thus having a 2^N complexity
            std::sort( next.begin(), next.end() );
            next.erase( std::unique( next.begin(), next.end() ), next.end() );
            current = next;
            next.clear();
        } while (current.size() > 0);

        // We remove duplicates and keep m_active_nodes sorted
        std::sort( m_active_nodes.begin(), m_active_nodes.end() );
        m_active_nodes.erase( std::unique( m_active_nodes.begin(), m_active_nodes.end() ), m_active_nodes.end() );

        // Then the active genes
        m_active_genes.clear();
        for (auto i = 0u; i<m_active_nodes.size(); ++i) 
        {
            if (m_active_nodes[i] >= m_n) 
            {
                unsigned int idx = (m_active_nodes[i] - m_n) * 3;
                m_active_genes.push_back(idx);
                m_active_genes.push_back(idx + 1);
                m_active_genes.push_back(idx + 2);
            }
        }
        for (auto i = 0u; i<m_m; ++i) 
        {
            m_active_genes.push_back(m_r * m_c * 3 + i);
        }
    }

private:
    // number of inputs
    unsigned int m_n;
    // number of outputs
    unsigned int m_m;
    // number of rows
    unsigned int m_r;
    // number of columns
    unsigned int m_c;
    // number of levels_back allowed
    unsigned int m_l;
    // the functions allowed
    std::vector<basis_function> m_f;
    // lower and upper bounds on all genes
    std::vector<unsigned int> m_lb;
    std::vector<unsigned int> m_ub;
    // active nodes idx (guaranteed to be always sorted)
    std::vector<unsigned int> m_active_nodes;
    // active genes idx
    std::vector<unsigned int> m_active_genes;
    // the actual expression encoded in a chromosome
    std::vector<unsigned int> m_x;
    // the random engine for the class
    std::default_random_engine m_e;
};



} // end of namespace dcgp

#endif // DCGP_EXPRESSION_H
