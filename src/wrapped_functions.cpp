#include <cmath>
#include <string>
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>

#include "wrapped_functions.h"
#include "std_overloads.h"


namespace dcgp {

double my_sum(double x, double y)
{
        return x + y;
}

double d_my_sum(const std::vector<double>& x, const std::vector<double>& y)
{
    unsigned int n = x.size() - 1;
    return x[n] + y[n];
}

std::string print_my_sum(const std::string& s1, const std::string& s2)
{
    if (s1 == s2) 
    {
        return "(2*"+s1+")";
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

double d_my_diff(const std::vector<double>& x, const std::vector<double>& y)
{
    unsigned int n = x.size() - 1;
    return x[n] - y[n];
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

double d_my_mul(const std::vector<double>& b, const std::vector<double>& c)
{
    unsigned int n = b.size() - 1u;
    double retval = 0.;
    for (auto j = 0u; j <= n; ++j) 
    {
        retval += b[n-j]*c[j];
    }
    return retval;
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
    else if (s1 == "1")
    {
        return s2;
    }
    else if (s2 == "1")
    {
        return s1;
    }
    return ("(" + s1 + "*" + s2 + ")");
}

double my_div(double x, double y)
{
        return x / y;
}

double d_my_div(const std::vector<double>& b, const std::vector<double>& c)
{
    unsigned int n = b.size() - 1u;
    std::vector<double> a(b.size());
    a[0] = b[0] / c[0];

    for (auto i = 1u; i <= n; ++i) 
    {
        a[i] = b[i];
        for (auto j = 1u; j <= i; ++j) 
        {
            a[i] -= c[j]*a[i-j];
        }
        a[i] /= c[0];
    }
    return a[n];
}

std::string print_my_div(const std::string& s1, const std::string& s2)
{
    if (s1 == "0" && s2 != "0")
    {
        return "0";
    }
    else if (s1 == s2)
    {
        return "1";
    }
    return ("(" + s1 + "/" + s2 + ")");
}

double my_pow(double b, double c)
{
        return pow(std::abs(b),c);
}

double d_my_pow(const std::vector<double>& b, const std::vector<double>& c)
{
    // We derive this by setting a = exp(c * ln(|b|))
    unsigned int n = b.size() - 1u;
    std::vector<double> a(b.size());
    std::vector<double> f(b.size());
    std::vector<double> g(b.size());

    // We take care of the abs
    double sign = 1;
    if (b[0] < 0)
    {
        sign = -1;
    }

    a[0] = pow(sign * b[0], c[0]);
    f[0] = c[0] * log(sign * b[0]);
    g[0] = log(sign * b[0]);

    // We start with g = log(b)
    for (auto i = 1u; i <= n; ++i)
    {
        g[i] = 0;
        for (auto j = 1u; j <= i-1; ++j)
        {
            g[i] += (i - j) * sign * b[j] * g[i - j];
        }    
        g[i] /= i;
        g[i] = (sign * b[i] - g[i]) / sign / b[0];
    }
    // Then we do f = c * g
    for (auto i = 1u; i <= n; ++i)
    {
        f[i] = 0;
        for (auto j = 0u; j <= i; ++j) 
        {
            f[i] += c[i-j] * g[j];
        }
    }
    // And finally a = exp(f)
    for (auto i = 1u; i <= n; ++i)
    {
        a[i] = 0;
        for (auto j = 0u; j <= i-1; ++j)
        {
            a[i] += (i - j) * a[j] * f[i - j];
        }
        a[i] /= i;
    }

    return a[n];
}

std::string print_my_pow(const std::string& s1, const std::string& s2)
{
    if (s1 == "0" && s2 != "0")
    {
        return "0";
    }
    else if (s1 == "1")
    {
        return "1";
    }
    else if (s2 == "0" && s1 != "0")
    {
        return "1";
    }
    else if (s2 == "1")
    {
        return s1;
    }
    return ("abs(" + s1 + ")^(" + s2 + ")");
}

double my_sqrt(double b, double c)
{
        return sqrt(std::abs(b));
}

double d_my_sqrt(const std::vector<double>& b, const std::vector<double>& c)
{
    (void)c;
    unsigned int n = b.size() - 1u;
    std::vector<double> a(b.size());
    double sign = 1;
    if (b[0] < 0)
    {
        sign = -1;
    }

    a[0] = sqrt(sign*b[0]);
    double alpha = 0.5;
    for (auto i = 1u; i <= n; ++i)
    {
        a[i] = 0;
        for (auto j = 0u; j <= i-1; ++j)
        {
            a[i] += (i * alpha - j * (alpha + 1)) * sign * b[i-j] * a[j];
        }    
        a[i] /= (i * sign * b[0]);
    }
    return a[n];
}

std::string print_my_sqrt(const std::string& s1, const std::string& s2)
{
    if (s1 == "0")
    {
        return "0";
    }
    else if (boost::algorithm::ends_with(s1, "^2"))
    {
        return "|" + s1 + "|";
    }
    return ("sqrt(|" + s1 + "|)");
}


} // dcgp namespace ends
