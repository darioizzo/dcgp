// reading a text file
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "../include/expression.hpp"
#include "../include/fitness_functions.hpp"
#include "../include/kernel_set.hpp"
#include "detail/es.hpp"
#include "detail/read_data.hpp"

using namespace dcgp;

int main() {
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
  es_params params{4, "active", 2, 0, 1000};
  es(in, out, ex, params);

  // We print out the final expression
  audi::stream(std::cout, "Final expression: ", ex(in_sym), "\n");
  audi::stream(std::cout, "Final value: ", quadratic_error(ex, in, out), "\n");
  return false;
}
