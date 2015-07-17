#include <iostream>
#include <iomanip>
#include <ctime>
#include "src/dcgp.h"


int main() {
    // Instantiate a cgp program with 1 inputs, 1 outputs, 1 row, 100 columns and 101 level-backs using the function set +,-,x,/ 
    // (1e-9 is the tolerance in case a hit count fitness is used, while 34534 is the seed)
    dcgp::program simple(3,4,2,3,4,dcgp::function_set::minimal,1e-9, 34534);
    std::cout << simple << std::endl;

    std::vector<double> in_num({2.,3.,4.});
    std::vector<std::string> in_sym({"x","y","z"});
    std::cout << "Numerical value = " << simple.compute_f(in_num) << std::endl;
    std::cout << "Symbolic value = " << simple.compute_f(in_sym) << std::endl;

    return 0;
}

/* Possible output:
CGP Encoding:
    Number of inputs:       3
    Number of outputs:      4
    Number of rows:         2
    Number of columns:      3
    Number of levels-back allowed:  4
    Tolerance (hit based fitness):  1e-09

    Resulting lower bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  ... ]
    Resulting upper bounds: [3, 2, 2, 3, 2, 2, 3, 4, 4, 3, 4, 4, 3, 6, 6, 3, 6, 6, 8, 8,  ... ]

    Current program (encoded):  [1, 1, 0, 0, 1, 1, 2, 4, 4, 2, 4, 2, 2, 6, 1, 2, 3, 5, 7, 2,  ... ]
    Active nodes:           [0, 1, 2, 4, 6, 7]
    Active genes:           [3, 4, 5, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21]

Numerical value = [72, 4, 6, 2]
Symbolic value = [(((2*y)*z)*y), z, (2*y), x]


**/
