#include <iostream>
#include <iomanip>
#include <ctime>
#include "src/dcgp.h"


int main() {
    // Instantiate a cgp expression with 1 inputs, 1 outputs, 1 row, 100 columns and 101 level-backs using the function set +,-,x,/ 
    // (1e-9 is the tolerance in case a hit count fitness is used, while 34534 is the seed)
    dcgp::expression simple(3,2,2,7,4,dcgp::function_set::minimal, 34534);
    std::cout << simple << std::endl;

    std::vector<double> in_num({2.,3.,4.});
    std::vector<std::string> in_sym({"x","y","z"});
    std::cout << "Point is:" << in_num << std::endl;
    std::cout << "Numerical value = " << simple(in_num) << std::endl;

    std::vector<std::vector<double> > jet_0 = simple.compute_d2(0,2,in_num);
    std::cout << "Numerical values d^n/dx^n = " << jet_0 << std::endl;

    std::vector<std::vector<double> > jet_1 = simple.compute_d2(1,2,in_num);
    std::cout << "Numerical values d^n/dy^n = " << jet_1 << std::endl;

    std::vector<std::vector<double> > jet_2 = simple.compute_d2(2,2,in_num);
    std::cout << "Numerical values d^n/dz^n = " << jet_2 << std::endl;

    std::cout << "Symbolic value = " << simple(in_sym) << std::endl;

    return 0;
}

/* Possible output:

d-CGP Expression:
    Number of inputs:       3
    Number of outputs:      2
    Number of rows:         2
    Number of columns:      7
    Number of levels-back allowed:  4

    Resulting lower bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  ... ]
    Resulting upper bounds: [3, 2, 2, 3, 2, 2, 3, 4, 4, 3, 4, 4, 3, 6, 6, 3, 6, 6, 3, 8,  ... ]

    Current expression (encoded):   [1, 2, 0, 3, 1, 1, 0, 1, 3, 2, 0, 3, 2, 6, 2, 0, 3, 4, 1, 6,  ... ]
    Active nodes:           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13]
    Active genes:           [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  ... ]

Point is:[2, 3, 4]
Numerical value = [0.20000000000000001, 48]
Numerical values d^n/dx^n = [[0.20000000000000001, 48], [0.23999999999999999, -16], [-0.30399999999999999, -24]]
Numerical values d^n/dy^n = [[0.20000000000000001, 48], [-0.040000000000000001, 0], [0.016, 0]]
Numerical values d^n/dz^n = [[0.20000000000000001, 48], [0.16, 52], [-0.064000000000000001, 36]]
Symbolic value = [(((x*(z-x))-((z-x)+1))/(y+(z-x))), (((x*(z-x))*z)*((z-x)+1))]

**/
