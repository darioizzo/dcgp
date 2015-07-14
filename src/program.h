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
    program(unsigned int n, 
            unsigned int m, 
            unsigned int r, 
            unsigned int c, 
            unsigned int l, 
            std::vector<basis_function> f, 
            unsigned int seed = rng::get_seed());

    void set(const std::vector<unsigned int> &x);
    
    template <class T>
    std::vector<T> compute_f(const std::vector<T>& in) const
    {  
        if(in.size() != m_n)
        {
            throw input_error("Input size is incompatible");
        }
        std::vector<unsigned int> to_evaluate(nodes_to_evaluate(m_x));
//for (auto i : to_evaluate) std::cout << " " << i; std::cout << std::endl;
        std::vector<T> retval(m_m);
        std::map<unsigned int, T> node;
        for (auto i : to_evaluate) {
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
    std::vector<unsigned int> nodes_to_evaluate(const std::vector<unsigned int>& x) const;

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
    // the actual program encoded in a chromosome
    std::vector<unsigned int> m_x;
    // the random engine for the class
    std::default_random_engine m_e;
};

std::ostream &operator<<(std::ostream &, const program &);

} // end of namespace dcgp

#endif // DCGP_PROGRAM_H
