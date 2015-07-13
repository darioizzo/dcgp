#include <iostream>
#include <iomanip>
#include "src/encoding.h"
#include "src/basis_function.h"
#include "src/wrapped_functions.h"
#include "src/std_overloads.h"


int main() {
    std::vector<dcgp::basis_function> f;
    f.emplace_back(dcgp::my_sum,dcgp::d_my_sum,dcgp::print_my_sum);
    f.emplace_back(dcgp::my_diff,dcgp::d_my_diff,dcgp::print_my_diff);
    f.emplace_back(dcgp::my_mul,dcgp::d_my_mul,dcgp::print_my_mul);
    f.emplace_back(dcgp::my_div,dcgp::d_my_div,dcgp::print_my_div);

    
    dcgp::encoding simple(2,4,2,3,4,f);

    std::vector<unsigned int> x({1,0,1,2,0,0,2,3,1,3,2,1,0,3,4,2,4,4,6,5,6,4});
    simple.is_valid(x);
    std::cout << simple.nodes_to_evaluate(x) << std::endl;
    std::cout << simple << std::endl;
    std::cout << simple.compute_f({2.,3.}, x) << std::endl;

    std::cout << simple.pretty({"x0","x1"}, x) << std::endl;

    return 0;
}
