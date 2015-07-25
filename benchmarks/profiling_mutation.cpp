#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include "../src/dcgp.h"


void perform_active_mutations(unsigned int in,
                  unsigned int out,
                  unsigned int rows,
                  unsigned int columns,
                  unsigned int levels_back,
                  unsigned int number_of_mutations,
                  std::vector<dcgp::basis_function> function_set)
{
    // Instatiate the expression
    dcgp::expression ex(in, out, rows, columns, levels_back, function_set, 123);
    clock_t begin = clock();
    for (auto i = 0u; i < number_of_mutations; ++i){
        ex.mutate_active();
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "In: " << in << ", Out: " << out << ", Rows: " << rows << ", Cols: " << columns << ", Levels-back: " << levels_back << std::endl;
    std::cout << number_of_mutations << " active mutations took: " << elapsed_secs << " seconds" << std::endl;
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
int main() {
    dcgp::function_set basic_set({"sum","diff","mul","div"});
    perform_active_mutations(2,4,2,3,4, 100000, basic_set());
    perform_active_mutations(2,4,10,10,11, 100000, basic_set());
    perform_active_mutations(2,4,20,20,21, 100000, basic_set());
    perform_active_mutations(1,1,1,100,101, 100000, basic_set());
    perform_active_mutations(1,1,2,100,101, 100000, basic_set());
    perform_active_mutations(1,1,3,100,101, 100000, basic_set());
    return 0;
}

