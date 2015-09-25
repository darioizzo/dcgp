#ifndef DCGP_EXPRESSION_H
#define DCGP_EXPRESSION_H

#include <vector>
#include <string>
#include <map>
#include <random>
#include <initializer_list>
#include <stdexcept>

#include "basis_function.h"
#include "rng.h"


namespace dcgp {

/// A d-CGP expression
/**
 * This class represent a mathematical expression as encoded using CGP and contains
 * algorithms that compute its value (numerical and symbolical) and its derivatives 
 * its fitness on a given input target set, as well as mutate the expression. 
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class expression {
public:             
    expression(unsigned int n, 
            unsigned int m, 
            unsigned int r, 
            unsigned int c, 
            unsigned int l, 
            std::vector<basis_function> f, 
            unsigned int seed = rng::get_seed()
            );

    void set(const std::vector<unsigned int> &x);
    void set(const std::initializer_list<unsigned int>& x)
    {
        set(std::vector<unsigned int>(x));
    }

    /// Gets the chromosome
    /** 
     * Gets the chromosome encoding the current expression
     *
     * \return The chromosome
    */
    const std::vector<unsigned int> & get() const {return m_x;};

    /// Gets the active genes
    /** 
     * Gets the idx of the active genes in the current chromosome (numbering is from 0)
     *
     * \return An std::vector containing the idx of the active genes in the current chromosome
    */
    const std::vector<unsigned int> & get_active_genes() const {return m_active_genes;};
    /// Gets the active nodes
    /** 
     * Gets the idx of the active nodes in the current chromosome.
     * The numbering starts from 0 at the first input node to then follow PPSN tutorial from Miller
     *
     * \return An std::vector containing the idx of the active nodes
    */
    const std::vector<unsigned int> & get_active_nodes() const {return m_active_nodes;};
    /// Gets the number of inputs
    /** 
     * Gets the number of inputs of the c_CGP expression
     *
     * \return the number of inputs
    */
    unsigned int get_n() const {return m_n;};
    /// Gets the number of outputs
    /** 
     * Gets the number of outputs of the c_CGP expression
     *
     * \return the number of outputs
    */
    unsigned int get_m() const {return m_m;};

    /// Gets the functions
    /** 
     * Gets the functions
     *
     * \return an std::vector<basis_function>
    */
    const std::vector<basis_function>& get_f() const {return m_f;};

    void mutate_active();
    
    template <class T>
    std::vector<T> operator()(const std::vector<T>& in) const
    {  
        if(in.size() != m_n)
        {
            throw std::invalid_argument("Input size is incompatible");
        }
//for (auto i : m_active_nodes) std::cout << " " << i; std::cout << std::endl;
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
//std::cout << i << ", " << node[i] << std::endl;
        }
        for (auto i = 0u; i<m_m; ++i)
        {
            retval[i] = node[m_x[(m_r * m_c) * 3 + i]];
        }
        return retval;
    }

    template <class T>
    std::vector<T> operator()(const std::initializer_list<T>& in) const
    {  
        std::vector<T> dummy(in);
        return (*this)(dummy);
    }

    std::vector<std::vector<double> > differentiate(unsigned int wrt, unsigned int degree, const std::vector<double>& in) const;
    std::string human_readable() const;

protected: 
    bool is_valid(const std::vector<unsigned int>& x) const;
    void update_active();

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

std::ostream &operator<<(std::ostream &, const expression &);

} // end of namespace dcgp

#endif // DCGP_EXPRESSION_H
