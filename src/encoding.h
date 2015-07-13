#ifndef DCGP_DCGP_H
#define DCGP_DCGP_H

#include "basis_function.h"
#include <vector>
#include <string>

namespace dcgp {

class encoding {
public:
    encoding(unsigned int n, unsigned int m, unsigned int c, unsigned int r, unsigned int l, std::vector<basis_function> f);
    bool is_valid(const std::vector<unsigned int>& x);
    std::vector<double> compute_f(const std::vector<double>& in, const std::vector<unsigned int>& x);
    std::string human_readable() const;

public: //TODO change to protected
    std::vector<unsigned int> nodes_to_evaluate(const std::vector<unsigned int>& x);

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
    // lower and upper bounds on all genes
    std::vector<unsigned int> m_lb;
    std::vector<unsigned int> m_ub;
};

std::ostream &operator<<(std::ostream &, const encoding &);

} // end of namespace dcgp

#endif // DCGP_DCGP_H
