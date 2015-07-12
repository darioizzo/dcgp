
#include <iostream>

#include "dcgp.h"


namespace dcgp {

    dcgp::dcgp(unsigned int n, unsigned int m, unsigned int c, unsigned int r, unsigned int l, std::vector<basis_function> f)
    {
        std::cout << "A" << f[0](1,2) << std::endl;
        std::cout << "B" << f[0].m_f(1,2) << std::endl;
        std::cout << "C" << f[0].m_df(0,1,2) << std::endl;
    }

} // end of namespace dcgp

