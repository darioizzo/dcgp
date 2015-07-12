#include <iostream>
#include <iomanip>
#include "src/dcgp.h"
#include "src/basis_function.h"

using namespace dcgp;

int main() {
    std::vector<basis_function> f;
    f.emplace_back(my_sum,d_my_sum);
    f.emplace_back(my_sum,d_my_sum);
    std::cout << f[0](1,2) << '\n';
    std::cout << f[0].m_df(0,1,2) << '\n';
    std::cout << f[0].m_df(1,1,2) << '\n';
    std::cout << f[0].m_df(2,1,2) << '\n';
    return 0;
}
