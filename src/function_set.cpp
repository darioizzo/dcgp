#include "function_set.h"
#include "wrapped_functions.h"

namespace dcgp {

std::vector<basis_function> minimal_set() 
{
    std::vector<basis_function> retval;
    retval.emplace_back(my_sum,d_my_sum,d_my_sum2,print_my_sum);
    retval.emplace_back(my_diff,d_my_diff,d_my_diff2,print_my_diff);
    retval.emplace_back(my_mul,d_my_mul,d_my_mul2,print_my_mul);
    retval.emplace_back(my_div,d_my_div,d_my_div2,print_my_div);
    return retval;
}
std::vector<basis_function> function_set::minimal = minimal_set();

} // end of namespace dcgp
