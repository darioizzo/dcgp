#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include <string>
#include <vector>
#include "exceptions.h"

namespace dcgp {

/*--------------------------------------------------------------------------
*                                  BINARY FUNCTIONS
*------------------------------------------------------------------------**/
// f = b + c
double my_sum(double b, double c);
double d_my_sum(const std::vector<double>& b, const std::vector<double>& c);
std::string print_my_sum(const std::string& s1, const std::string& s2);

// f = b - c
double my_diff(double b, double c);
double d_my_diff(const std::vector<double>& b, const std::vector<double>& c);
std::string print_my_diff(const std::string& s1, const std::string& s2);

// f = b * c
double my_mul(double b, double c);
double d_my_mul(const std::vector<double>& b, const std::vector<double>& c);
std::string print_my_mul(const std::string& s1, const std::string& s2);

// f = b / c
double my_div(double b, double c);
double d_my_div(const std::vector<double>& b, const std::vector<double>& c);
std::string print_my_div(const std::string& s1, const std::string& s2);

// f = pow(|b|,c)
double my_pow(double b, double c);
double d_my_pow(const std::vector<double>& b, const std::vector<double>& c);
std::string print_my_pow(const std::string& s1, const std::string& s2);

/*--------------------------------------------------------------------------
*                                  UNARY FUNCTIONS
*------------------------------------------------------------------------**/

// f = sqrt(|b|)
double my_sqrt(double b, double c);
double d_my_sqrt(const std::vector<double>& b, const std::vector<double>& c);
std::string print_my_sqrt(const std::string& s1, const std::string& s2);

/*--------------------------------------------------------------------------
*                                  HELPER FUNCTIONS
*------------------------------------------------------------------------**/

double d_not_implemented(const std::vector<double>& b, const std::vector<double>& c);

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H