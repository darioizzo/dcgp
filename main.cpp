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

/* Possible output:
CGP Encoding:
    Number of inputs:       3
    Number of outputs:      1
    Number of rows:         5
    Number of columns:      5
    Number of levels-back allowed:  6

    Resulting lower bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  ... ]
    Resulting upper bounds: [3, 2, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2, 3, 7, 7, 3, 7, 7, 3, 7, 7, 3, 7, 7, 3, 7, 7, 3, 12, 12, 3, 12, 12, 3, 12, 12, 3, 12, 12, 3, 12, 12, 3, 17, 17, 3, 17, 17, 3, 17, 17, 3, 17, 17, 3, 17, 17,  ... ]

    Current program (encoded):  [1, 2, 0, 2, 2, 1, 0, 1, 2, 3, 2, 2, 2, 0, 1, 2, 6, 7, 2, 4, 7, 3, 0, 5, 0, 1, 6, 1, 0, 6, 2, 11, 3, 0, 10, 11, 3, 6, 9, 1, 8, 5, 2, 2, 9, 2, 0, 11, 0, 11, 10, 3, 11, 13, 2, 3, 6, 0, 10, 17,  ... ]

Numerical value = [0.5]
Symbolic value = [((y+1)/((y+1)*(z-x)))]

**/
