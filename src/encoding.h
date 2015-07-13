#ifndef DCGP_DCGP_H
#define DCGP_DCGP_H

#include <vector>
#include <string>
#include <map>

#include "basis_function.h"
#include "exceptions.h"


namespace dcgp {

class encoding {
public:
    encoding(unsigned int n, unsigned int m, unsigned int r, unsigned int c, unsigned int l, std::vector<basis_function> f);
    bool is_valid(const std::vector<unsigned int>& x) const;
    
    template <class T>
    std::vector<T> compute_f(const std::vector<T>& in, const std::vector<unsigned int>& x) const
    {  
        if (!is_valid(x)) throw input_error("Invalid chromosome");
        std::vector<unsigned int> to_evaluate(nodes_to_evaluate(x));
        std::vector<T> retval(m_m);
        std::map<unsigned int, T> node;
        for (auto i : to_evaluate) {
            if (i < m_n) 
            {
                node[i] = in[i];
            } else {
                unsigned int idx = (i - 2) * 3;
                node[i] = m_f[x[idx]](node[x[idx + 1]], node[x[idx + 2]]);
            }
        }
        for (auto i = 0u; i<m_m; ++i)
        {
            retval[i] = node[x[(m_r * m_c) * 3 + i]];
        }
        return retval;
    }
    //std::vector<std::string>  pretty(const std::vector<std::string>& in, const std::vector<unsigned int>& x) const;
    std::string human_readable() const;

private: 
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
};

std::ostream &operator<<(std::ostream &, const encoding &);

} // end of namespace dcgp

#endif // DCGP_DCGP_H
