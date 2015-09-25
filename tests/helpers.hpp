#ifndef DCGP_HELPERS_FUNCTION_H
#define DCGP_HELPERS_FUNCTION_H

#include <vector>
#include <stdexcept>

namespace dcgp
{

template <typename Iterator>
inline bool EPSILON_COMPARE(Iterator It1first, Iterator It1last, Iterator It2first, double epsilon)
{
    return std::equal(It1first, It1last, It2first, [epsilon](double x, double y){return std::abs(x-y) < epsilon;});
}

template <typename T>
inline bool EPSILON_COMPARE(std::vector<T> first, std::vector<T> second, double epsilon)
{
    if (first.size()!=second.size()) {
        throw std::invalid_argument("Vectors must be of equal size");
    }

    return EPSILON_COMPARE(first.begin(), first.end(), second.begin(), epsilon);
}



} // end of namespace dcgp 

#endif
