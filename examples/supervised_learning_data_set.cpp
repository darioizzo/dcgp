// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include <sstream>

#include "../include/function_set.hpp"
#include "../include/expression.hpp"
#include "../include/fitness_functions.hpp"
#include "../include/io.hpp"
#include "detail/read_data.hpp"
#include "detail/es.hpp"

using namespace dcgp;

int main () {
    // Random seed
    std::random_device rd;
    // We read the data from file
    std::vector<std::vector<double> > in, out;
    read_data(in, out, "../../examples/data/pressure.data");

    // Function set
    dcgp::function_set basic_set({"sum", "diff", "mul", "div", "sig"});

    // d-CGP expression
    dcgp::expression ex(2, 1, 2, 100, 101, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"x","y"});

    // We use a simple ES(1+4) to evolve an expression that represents our target.
    // Mutation only mutates 2 active genes
    es_params params{4, "active", 4, 0.01, 100000};
    es(in, out, ex, params);

    // We print out the final expression
    stream(std::cout, "Final expression: ", ex(in_sym), "\n");
    stream(std::cout, "Final value: ", quadratic_error(ex, in, out), "\n");
    return false;

}
