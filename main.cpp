#include <iostream>
#include <iomanip>
#include <ctime>
#include "src/dcgp.h"


int main() {
    // Instantiate a cgp expression with 1 inputs, 1 outputs, 1 row, 100 columns and 101 level-backs using the function set +,-,x,/ 
    // (1e-9 is the tolerance in case a hit count fitness is used, while 34534 is the seed)
    dcgp::expression simple(3,4,2,3,4,dcgp::function_set::minimal,1e-9, 34534);
    std::cout << simple << std::endl;

    std::vector<double> in_num({2.,3.,4.});
    std::vector<std::string> in_sym({"x","y","z"});
    std::cout << "Point is:" << in_num << std::endl;
    std::cout << "Numerical value = " << simple.compute(in_num) << std::endl;
    std::cout << "Numerical value d/dx = " << simple.compute_d(0,in_num) << std::endl;
    std::cout << "Numerical value d/dy = " << simple.compute_d(1,in_num) << std::endl;
    std::cout << "Numerical value d/dz = " << simple.compute_d(2,in_num) << std::endl;
    std::cout << "Symbolic value = " << simple.compute(in_sym) << std::endl;

    return 0;
}

/* Possible output:

d-CGP Expression:
    Number of inputs:       3
    Number of outputs:      4
    Number of rows:         2
    Number of columns:      3
    Number of levels-back allowed:  4
    Tolerance (hit based fitness):  1e-09

    Resulting lower bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  ... ]
    Resulting upper bounds: [3, 2, 2, 3, 2, 2, 3, 4, 4, 3, 4, 4, 3, 6, 6, 3, 6, 6, 8, 8,  ... ]

    Current expression (encoded):   [1, 2, 0, 3, 1, 1, 0, 1, 3, 2, 0, 3, 2, 6, 2, 0, 3, 4, 5, 6,  ... ]
    Active nodes:           [0, 1, 2, 3, 4, 5, 6, 8]
    Active genes:           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18, 19, 20, 21]

Numerical value = [5, 4, 3, 4]
Symbolic value = [(y+(z-x)), (x*(z-x)), ((z-x)+1), (x*(z-x))]

**/
