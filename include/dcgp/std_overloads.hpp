#ifndef STD_OVERLOADS_H
#define STD_OVERLOADS_H

#include <vector>

template <typename T>
inline std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
    os << "[";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {
        os << " " << *ii;
    }
    os << " ]";
    return os;
}

#endif
