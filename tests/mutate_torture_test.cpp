#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include "../src/dcgp.h"


bool mutate_a_lot(unsigned int in,
                  unsigned int out,
                  unsigned int rows,
                  unsigned int columns,
                  unsigned int levels_back,
                  std::vector<dcgp::basis_function> function_set)
{
    // Instatiate the expression
    dcgp::expression ex(in, out, rows, columns, levels_back, function_set);
    // Pick a radnom input value in [-1, 1)
    std::vector<double> in_num;
    std::default_random_engine re;
    for (auto i = 0u; i < in; ++i)
    {
        in_num.push_back(std::uniform_real_distribution<double>(-1, 1)(re));
    }
    clock_t begin = clock();
    double fit;
    for (auto i = 0u; i < 10000; ++i){
        ex.mutate();
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "In: " << in << ", Out: " << out << ", Rows: " << rows << ", Cols: " << columns << ", Levels-back: " << levels_back << std::endl;
    std::cout << "10000 mutations took: " << elapsed_secs << " seconds" << std::endl;
    return false;
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
int main() {
    std::vector<dcgp::basis_function> function_set = dcgp::function_set::minimal;
    return mutate_a_lot(2,4,2,3,4, function_set) ||
           mutate_a_lot(2,4,10,10,11, function_set) ||
           mutate_a_lot(2,4,20,20,21, function_set) ||
           mutate_a_lot(1,1,1,100,101, function_set) ||
           mutate_a_lot(1,1,2,100,101, function_set) ||
           mutate_a_lot(1,1,3,100,101, function_set);
}

