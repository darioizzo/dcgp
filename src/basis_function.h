#ifndef DCGP_BASIS_FUNCTION_H
#define DCGP_BASIS_FUNCTION_H

namespace dcgp {

class basis_function {
public:
    virtual double operator()(double a, double b) = 0;
    virtual double d1(double a, double b) = 0;
    virtual double d2(double a, double b) = 0;
};

} // end of namespace dcgp

#endif // DCGP_BASIS_FUNCTION_H
