#include <iostream>
#include <iomanip>
#include "src/program.h"
#include "src/basis_function.h"
#include "src/wrapped_functions.h"
#include "src/std_overloads.h"


int main() {
    std::vector<dcgp::basis_function> f;
    f.emplace_back(dcgp::my_sum,dcgp::d_my_sum,dcgp::print_my_sum);
    f.emplace_back(dcgp::my_diff,dcgp::d_my_diff,dcgp::print_my_diff);
    f.emplace_back(dcgp::my_mul,dcgp::d_my_mul,dcgp::print_my_mul);
    f.emplace_back(dcgp::my_div,dcgp::d_my_div,dcgp::print_my_div);

    dcgp::program simple(4,1,3,6,3,f,9);

    //std::vector<unsigned int> x({0, 0, 1, 1, 0, 0, 1, 3, 1, 2, 0, 1, 0, 4, 4, 2, 5, 4, 2, 5, 7, 3});
    std::cout << simple << std::endl;
    std::cout << simple.compute_f(std::vector<double>({2.,3.,4.,-2.})) << std::endl;
    std::cout << simple.compute_f(std::vector<std::string>({"x","y","z","t"})) << std::endl;
    //simple.set(x);
    //std::cout << simple << std::endl;
    //std::cout << simple.compute_f(std::vector<double>({2.,3.})) << std::endl;
    //std::cout << simple.compute_f(std::vector<std::string>({"x0","x1"})) << std::endl;


    return 0;
}

