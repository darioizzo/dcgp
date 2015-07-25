#include <string>
#include <iostream>

#include "basis_function.h"

namespace dcgp {

std::ostream& operator<<(std::ostream& os, const basis_function& obj)
{
    os << obj.m_name;
    return os;
}

} // end of namespace dcgp

 