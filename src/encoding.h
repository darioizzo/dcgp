#ifndef DCGP_DCGP_H
#define DCGP_DCGP_H

#include "basis_function.h"
#include <vector>

namespace dcgp {

class encoding {
public:
    encoding(unsigned int n, unsigned int m, unsigned int c, unsigned int r, unsigned int l, std::vector<basis_function> f);
    

    bool is_valid(std::vector<unsigned int> x, bool verbose = false);
private:
    // number of inputs
    unsigned int m_n;
    // number of outputs
    unsigned int m_m;
    // number of columns
    unsigned int m_c;
    // number of rows
    unsigned int m_r;
    // number of levels_back allowed
    unsigned int m_l;
    // the functions allowed
    std::vector<basis_function> m_f;
    std::vector<unsigned int> m_lb;
    std::vector<unsigned int> m_ub;
};

} // end of namespace dcgp

#endif // DCGP_DCGP_H
