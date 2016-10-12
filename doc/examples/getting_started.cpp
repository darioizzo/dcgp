#include <audi/audi.hpp>
#include <iostream>

using namespace audi;

int main() {
    // We want to compute the Taylor expansion of a function f (and thus all derivatives) at x=2, y=3
    using gdual = gdual<double>;
    // 1 - Define the generalized dual numbers (over doubles, 7 is the truncation order, i.e. the maximum
    // order of derivation we will need)
    gdual x(2, "x", 7);
    gdual y(3, "y", 7);

    // 2 - Compute your function as usual
    auto f = exp(x*x + cbrt(y) / log(x*y));

    // 3 - Inspect the results (this has a constant complexity now as all computations have been made already)
    std::cout << "Taylor polynomial: " << f << std::endl;                      // This is the Taylor expansion of f (truncated at the 7th order)
    std::cout << "Derivative value: " << f.get_derivative({1,0}) << std::endl; // This is the value of the derivative (d / dx)
    std::cout << "Derivative value: " << f.get_derivative({4,3}) << std::endl; // This is the value of the mixed derivative (d^7 / dx^4dy^3)

    // 4 - Using the dictionary interface (note the presence of the "d" before all variables)
    std::cout << "Derivative value: " << f.get_derivative({{"dx", 1}}) << std::endl; // This is the value of the derivative (d / dx)
    std::cout << "Derivative value: " << f.get_derivative({{"dx", 4}, {"dy", 3}}) << std::endl; // This is the value of the mixed derivative (d^7 / dx^4dy^3)
}
