#ifndef DCGP_FUNCTION_SET_H
#define DCGP_FUNCTION_SET_H

#include <vector>
#include "basis_function.h"

namespace dcgp {

class function_set
{
public:
    static std::vector<dcgp::basis_function> minimal;
};

} // end of namespace dcgp

#endif // DCGP_FUNCTION_SET_H
 