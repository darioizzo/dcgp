#include <audi/audi.hpp>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <random>

#include "include/expression.hpp"
#include "include/kernel.hpp"
#include "include/kernel_set.hpp"
#include "include/wrapped_functions.hpp"

using namespace dcgp;
using fun_type = std::function<double(const std::vector<double> &)>;
using fun_print_type =
    std::function<std::string(const std::vector<std::string> &)>;

int main() {
  expression<gdual_d> ex(1, 1, 1, 6, 6, 2,
                         kernel_set<gdual_d>({"sum", "diff", "div", "mul"})(),
                         67);
  gdual_d x(1.2, "x", 5);
  audi::stream(std::cout, ex({x}), "\n");
}
