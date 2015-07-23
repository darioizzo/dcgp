#ifndef DCGP_FUNCTION_SET_H
#define DCGP_FUNCTION_SET_H

#include <vector>
#include "basis_function.h"

namespace dcgp {

/// Function set
/**
 * Contains, as static members, several std::vector of dcgp::basis_function containing 
 * function sets of common use. The user can access each set via the syntax 
 * dcgp::function_set::SET_NAME
 *
 * @author Dario Izzo (dario.izzo@gmail.com)
 */
class function_set
{
public:
	/// The minimal function set, containing only +,-,*,/
    static std::vector<dcgp::basis_function> minimal;
    static std::vector<dcgp::basis_function> extended;
};

} // end of namespace dcgp

#endif // DCGP_FUNCTION_SET_H
 