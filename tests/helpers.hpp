#ifndef DCGP_HELPERS_FUNCTION_H
#define DCGP_HELPERS_FUNCTION_H

#include <boost/test/included/unit_test.hpp>
#include <stdexcept>
#include <vector>

namespace dcgp
{

template <typename T>
void CHECK_EQUAL_V(const T &in1, const T &in2)
{
    BOOST_CHECK_EQUAL(in1.size(), in2.size());
    for (decltype(in1.size()) i = 0u; i < in1.size(); ++i) {
        BOOST_CHECK_EQUAL(in1[i], in2[i]);
    }
}

template <typename T>
void CHECK_CLOSE_V(const T &in1, const T &in2, double tol)
{
    BOOST_CHECK_EQUAL(in1.size(), in2.size());
    for (decltype(in1.size()) i = 0u; i < in1.size(); ++i) {
        BOOST_CHECK_CLOSE(in1[i], in2[i], tol);
    }
}

} // end of namespace dcgp

#endif
