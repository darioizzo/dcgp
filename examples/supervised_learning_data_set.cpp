// reading a text file
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/fitness_functions.hpp>
#include <dcgp/kernel_set.hpp>

#include "detail/es.hpp"
#include "detail/read_data.hpp"

using namespace dcgp;

int main()
{
    // Random seed
    std::random_device rd;
    // We read the data from file
    std::vector<std::vector<double>> in, out;
    read_data(in, out, "../../examples/data/symbolic.data");

    // Function set
    dcgp::kernel_set<double> basic_set({"sum", "diff", "mul", "div"});

    // d-CGP expression
    dcgp::expression<double> ex(1, 1, 1, 15, 16, 2, basic_set(), rd());

    // Symbols
    std::vector<std::string> in_sym({"x"});

    // We use a simple ES(1+4) to evolve an expression that represents our target.
    // Mutation only mutates 2 active genes
    es_params params{4, "active", 2, 0, 100000};
    es(in, out, ex, params);

    // We print out the final expression
    audi::print("Final expression: ", ex(in_sym), "\n");
    audi::print("Final value: ", ex.loss(in, out, "MSE", true), "\n");
    return false;
}
