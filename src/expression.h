#ifndef DCGP_EXPRESSION_H
#define DCGP_EXPRESSION_H

#include <vector>
#include <string>
#include <map>
#include <random>

#include "basis_function.h"
#include "exceptions.h"
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
    const std::vector<unsigned int> & get() const;

    const std::vector<unsigned int> & get_active_genes() const;
    const std::vector<unsigned int> & get_active_nodes() const;
    unsigned int get_n() const;
    unsigned int get_m() const;

    void mutate();
    
    template <class T>
    std::vector<T> operator()(const std::vector<T>& in) const
    {  
        if(in.size() != m_n)
        {
            throw input_error("Input size is incompatible");
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

    std::vector<double> compute_d(unsigned int wrt, const std::vector<double>& in) const;
    std::vector<std::vector<double> > compute_d2(unsigned int wrt, unsigned int degree, const std::vector<double>& in) const;
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
