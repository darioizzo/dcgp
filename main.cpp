#include <iostream>
#include <iomanip>
#include "src/dcgp.h"


int main() {
    // Instantiate a cgp program with 3 inputs, 2 outputs, 1 row, 5 columns and 6 level-backs using the function set +,-,x,/
    dcgp::program simple(3,1,5,5,6,dcgp::function_set::minimal);
    std::cout << simple << std::endl;

    std::vector<double> in_num({2.,3,4});
    std::vector<std::string> in_sym({"x","y","z"});
    std::cout << "Numerical value = " << simple.compute_f(in_num) << std::endl;
    std::cout << "Symbolic value = " << simple.compute_f(in_sym) << std::endl;
    return 0;
}

