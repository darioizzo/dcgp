#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <string>
#include <vector>
#include "exceptions.h"

namespace dcgp {

double my_sum(double x, double y);
double d_my_sum(unsigned int i, double x, double y);
double d_my_sum2(const std::vector<double>& x, const std::vector<double>& y);
std::string print_my_sum(const std::string& s1, const std::string& s2);

double my_diff(double x, double y);
double d_my_diff(unsigned int i, double x, double y);
double d_my_diff2(const std::vector<double>& x, const std::vector<double>& y);
std::string print_my_diff(const std::string& s1, const std::string& s2);

double my_mul(double x, double y);
double d_my_mul(unsigned int i, double x, double y);
double d_my_mul2(const std::vector<double>& x, const std::vector<double>& y);
std::string print_my_mul(const std::string& s1, const std::string& s2);

double my_div(double x, double y);
double d_my_div(unsigned int i, double x, double y);
double d_my_div2(const std::vector<double>& x, const std::vector<double>& y);
std::string print_my_div(const std::string& s1, const std::string& s2);

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H