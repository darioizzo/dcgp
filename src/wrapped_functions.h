#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <string>
#include "exceptions.h"

namespace dcgp {

double my_sum(double x, double y)
{
        return x + y;
}

double d_my_sum(unsigned int i, double x, double y)
{
    (void)x; (void)y;
    if (i < 2) return 1;
    else throw derivative_error("Index is out of bounds");
}

std::string print_my_sum(const std::string& s1, const std::string& s2)
{
    if (s1 == s2) 
    {
        return "2"+s1;
    }
    else if (s1 == "0")
    {
        return s2;
    }
    else if (s2 == "0")
    {
        return s1;
    }
    return ("(" + s1 + "+" + s2 + ")");
}


double my_diff(double x, double y)
{
        return x - y;
}

double d_my_diff(unsigned int i, double x, double y)
{
    (void)x; (void)y;
    if (i == 0) return 1;
    else if (i==1) return -1;
    else throw derivative_error("Index is out of bounds");
}

std::string print_my_diff(const std::string& s1, const std::string& s2)
{
    if (s1 == s2) 
    {
        return "0";
    }
    else if (s1 == "0")
    {
        return "(-" + s2 + ")";
    }
    else if (s2 == "0")
    {
        return s1;
    }
    return ("(" + s1 + "-" + s2 + ")");
}

double my_mul(double x, double y)
{
        return (x * y);
}

double d_my_mul(unsigned int i, double x, double y)
{
    if (i == 0) return y;
    else if (i==1) return x;
    else throw derivative_error("Index is out of bounds");
}

std::string print_my_mul(const std::string& s1, const std::string& s2)
{
    if (s1 == "0" || s2 == "0")
    {
        return "0";
    }
    else if (s1 == s2)
    {
        return s1 + "^2";
    }
    return ("(" + s1 + "*" + s2 + ")");
}

double my_div(double x, double y)
{
        return x * y;
}

double d_my_div(unsigned int i, double x, double y)
{
    if (i == 0) return 1. / y;
    else if (i==1) return - x / y / y;
    else throw derivative_error("Index is out of bounds");
}

std::string print_my_div(const std::string& s1, const std::string& s2)
{
    return ("(" + s1 + "/" + s2 + ")");
}

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H