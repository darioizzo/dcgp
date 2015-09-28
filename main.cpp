#include <iostream>
#include <iomanip>
#include <ctime>
#include "src/dcgp.h"


int main() {
    // We define the set of functions we want to use
    dcgp::function_set basic_set({"sum","diff","mul","div"});

    // We instantiate a d-CGP expression
    unsigned int n_inputs = 3;
    unsigned int n_outputs = 1;
    unsigned int n_rows = 1;
    unsigned int n_columns = 50;
    unsigned int n_level_backs = 51;
    dcgp::expression simple(n_inputs,n_outputs,n_rows,n_columns,n_level_backs,basic_set());

    // We inspect it
    std::cout << simple << std::endl;

    // We compute the expression value and its derivatives in a point
    std::vector<double> in_num({2.,3.,4.});
    std::cout << "Point is:" << in_num << std::endl;
    std::cout << "Numerical value = " << simple(in_num) << std::endl;

    std::vector<audi::gdual> jet_0 = simple.differentiate(in_num,1);
    std::cout << "Numerical values d/dx = " << jet_0[0].get_derivative({1,0,0}) << std::endl;
    std::cout << "Numerical values d/dy = " << jet_0[0].get_derivative({0,1,0}) << std::endl;
    std::cout << "Numerical values d/dz = " << jet_0[0].get_derivative({0,0,1}) << std::endl;
    /*
    std::vector<std::vector<double> > jet_1 = simple.differentiate(1,2,in_num);
    std::cout << "Numerical values d^n/dy^n = " << jet_1 << std::endl;

    std::vector<std::vector<double> > jet_2 = simple.differentiate(2,2,in_num);
    std::cout << "Numerical values d^n/dz^n = " << jet_2 << std::endl;
    */

    // We stream a symbolic representation of the expression
    std::vector<std::string> in_sym({"x","y","z"});
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
