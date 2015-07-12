#ifndef DCGP_WRAPPED_FUNCTIONS_H
#define DCGP_WRAPPED_FUNCTIONS_H

#include "exceptions.h"

namespace dcgp {

double my_sum(double x, double y)
{
        return x + y;
}

double d_my_sum(unsigned int i, double x, double y)
{
	if (i < 2) return 1;
	else throw derivative_error("Index is out of bounds");
}

double my_diff(double x, double y)
{
        return x - y;
}

double d_my_diff(unsigned int i, double x, double y)
{
	if (i == 0) return 1;
	else if (i==1) return -1;
	else throw derivative_error("Index is out of bounds");
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

} // dcgp namespace ends

#endif // DCGP_WRAPPED_FUNCTIONS_H