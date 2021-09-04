#ifndef DCGP_HELPERS_FUNCTION_H
#define DCGP_HELPERS_FUNCTION_H

#include <boost/test/included/unit_test.hpp>
#include <stdexcept>
#include <vector>

namespace dcgp
{

// Have to make it a macro so that it reports exact line numbers when checks fail.
#define CHECK_EQUAL_V(in1, in2) { \
    BOOST_CHECK_EQUAL(in1.size(), in2.size()); \
    for (decltype(in1.size()) i = 0u; i < in1.size(); ++i) { \
        BOOST_CHECK_EQUAL(in1[i], in2[i]); \
    } \
}

#define CHECK_CLOSE_V(in1, in2, tol) { \
    BOOST_CHECK_EQUAL(in1.size(), in2.size()); \
    for (decltype(in1.size()) i = 0u; i < in1.size(); ++i) { \
        BOOST_CHECK_CLOSE(in1[i], in2[i], tol); \
    } \
}



} // end of namespace dcgp

#endif
