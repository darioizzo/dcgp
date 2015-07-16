#include <iostream>
#include <iomanip>
#include <ctime>
#include "src/dcgp.h"


int main() {
    // Instantiate a cgp program with 3 inputs, 2 outputs, 1 row, 5 columns and 6 level-backs using the function set +,-,x,/
    dcgp::program simple(1,1,1,100,101,dcgp::function_set::minimal,34534);
    std::cout << simple << std::endl;

    std::vector<double> in_num({2.});
    std::vector<std::string> in_sym({"x"});
    std::cout << "Numerical value = " << simple.compute_f(in_num) << std::endl;
    std::cout << "Symbolic value = " << simple.compute_f(in_sym) << std::endl;

    clock_t begin = clock();
    for (auto i = 0u; i < 10000; ++i){
        simple.mutate();
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "10000 mutations took: " << elapsed_secs << " seconds" << std::endl;

    return 0;
}

/* Possible output:
CGP Encoding:
    Number of inputs:       3
    Number of outputs:      4
    Number of rows:         1
    Number of columns:      5
    Number of levels-back allowed:  6

    Resulting lower bounds: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Resulting upper bounds: [3, 2, 2, 3, 3, 3, 3, 4, 4, 3, 5, 5, 3, 6, 6, 7, 7, 7, 7]

    Current program (encoded):  [1, 1, 2, 2, 0, 0, 0, 0, 2, 1, 5, 4, 0, 3, 5, 1, 6, 0, 0]

    Active nodes:   [0, 1, 2, 4, 5, 6]

    Active genes:   [3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18]

Numerical value = [3, 2, 2, 2]
Symbolic value = [y, ((x+t)-x^2), x, x]
1000000 mutations took: 0.014231 seconds
**/
