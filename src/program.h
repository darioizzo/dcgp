#ifndef DCGP_PROGRAM_H
#define DCGP_PROGRAM_H

#include <vector>
#include <string>
#include <map>
#include <random>

#include "basis_function.h"
#include "exceptions.h"
#include "rng.h"


namespace dcgp {

class program {
public:

    enum fitness_type { ERROR_BASED, HITS_BASED };

    program(unsigned int n, 
            unsigned int m, 
            unsigned int r, 
            unsigned int c, 
            unsigned int l, 
            std::vector<basis_function> f, 
            double tol = 1e-12,
            unsigned int seed = rng::get_seed()
            );

    void set(const std::vector<unsigned int> &x);
    const std::vector<unsigned int> & get() const;

    const std::vector<unsigned int> & get_active_genes() const;
    const std::vector<unsigned int> & get_active_nodes() const;

    void mutate();
    double fitness(const std::vector<std::vector<double> >& in_des, const std::vector<std::vector<double> >& out_des, fitness_type type = fitness_type::ERROR_BASED) const;
    
    template <class T>
    std::vector<T> compute_f(const std::vector<T>& in) const
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
    // active nodes
    std::vector<unsigned int> m_active_nodes;
    // active genes
    std::vector<unsigned int> m_active_genes;
    // the actual program encoded in a chromosome
    std::vector<unsigned int> m_x;
    // tolerance for the hits based fitness
    double m_tol;
    // the random engine for the class
    std::default_random_engine m_e;
};

std::ostream &operator<<(std::ostream &, const program &);

} // end of namespace dcgp

#endif // DCGP_PROGRAM_H
