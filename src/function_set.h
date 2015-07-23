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
	function_set();
	function_set(const std::vector<std::string>&);
	void push_back(const std::string&);
	std::vector<dcgp::basis_function> operator()() const;
private:
    std::vector<dcgp::basis_function> m_functions;
};

} // end of namespace dcgp

#endif // DCGP_FUNCTION_SET_H
 