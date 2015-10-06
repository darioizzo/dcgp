// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include <sstream>

#include "../src/function_set.hpp"
#include "../src/expression.hpp"
#include "../src/fitness_functions.hpp"
#include "../src/std_overloads.hpp"
#include "detail/read_data.hpp"
#include "detail/es.hpp"

using namespace std;

int main () {
    // Random seed
    std::random_device rd;
    // We read the data from file
    std::vector<std::vector<double> > in, out;
    read_data(in, out, "../../examples/data/symbolic.data");

    // Function set
    dcgp::function_set basic_set({"sum", "diff", "mul", "div", "sin"});

    // d-CGP expression
    dcgp::expression ex(1, 1, 1, 15, 16, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"x"});

    // We use a simple ES(1+4) to evolve an expression that represents our target. 
    // Mutation only mutates 2 active genes
    es_params params{4, "active", 2, 0};
    es(in, out, ex, params);

    // We print out the final expression
    std::cout << "Final expression: " << ex(in_sym) << std::endl;
    return false;
    
}

