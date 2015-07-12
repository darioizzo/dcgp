#ifndef DCGP_DCGP_H
#define DCGP_DCGP_H

#include "basis_function.h"
#include <vector>

namespace dcgp {

class dcgp {
public:
    dcgp(unsigned int n, unsigned int m, unsigned int c, unsigned int r, unsigned int l, std::vector<basis_function> f);
};

} // end of namespace dcgp

#endif // DCGP_DCGP_H
