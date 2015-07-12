#include <iostream>
#include <iomanip>
#include "src/encoding.h"
#include "src/basis_function.h"
#include "src/wrapped_functions.h"


int main() {
    std::vector<dcgp::basis_function> f;
    f.emplace_back(dcgp::my_sum,dcgp::d_my_sum);
    std::cout << f[0](1,2) << '\n';
    std::cout << f[0].m_df(0,1,2) << '\n';
    std::cout << f[0].m_df(1,1,2) << '\n';
    
    dcgp::encoding simple(2,2,2,2,1,f);

    simple.is_valid({0,1,1,0,1,1,0,3,3,0,3,3,4,5}, true);

    return 0;
}
