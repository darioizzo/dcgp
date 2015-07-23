#include "function_set.h"
#include "wrapped_functions.h"

namespace dcgp {

std::vector<basis_function> minimal_set() 
{
    std::vector<basis_function> retval;
    retval.emplace_back(my_sum,d_my_sum,print_my_sum);
    retval.emplace_back(my_diff,d_my_diff,print_my_diff);
    retval.emplace_back(my_mul,d_my_mul,print_my_mul);
    retval.emplace_back(my_div,d_my_div,print_my_div);
    return retval;
}
std::vector<basis_function> function_set::minimal = minimal_set();

std::vector<basis_function> extended_set() 
{
    std::vector<basis_function> retval;
    retval.emplace_back(my_sum,d_my_sum,print_my_sum);
    retval.emplace_back(my_diff,d_my_diff,print_my_diff);
    retval.emplace_back(my_mul,d_my_mul,print_my_mul);
    retval.emplace_back(my_div,d_my_div,print_my_div);
    retval.emplace_back(my_sqrt,d_my_sqrt,print_my_sqrt);
    return retval;
}
std::vector<basis_function> function_set::extended = extended_set();

} // end of namespace dcgp
