#include <string>
#include <iostream>

#include "basis_function.h"

namespace dcgp {


/// Overload stream operator for dcgp::basis_function
/**
 * @param[in] obj dcgp::basis_function to be inserted into the stream.
 * @param[out] os std::ostream to which the problem will be streamed.
 *
 * @return reference to os.
 */
std::ostream& operator<<(std::ostream& os, const basis_function& obj)
{
    os << obj.m_name;
    return os;
}

} // end of namespace dcgp

 