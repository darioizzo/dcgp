#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include <audi/audi.hpp>

#include "include/basis_function.hpp"
#include "include/expression.hpp"
#include "include/wrapped_functions.hpp"
#include "include/function_set.hpp"
#include "include/io.hpp"

using namespace dcgp;
using fun_type = std::function<double(const std::vector<double>&)>;
using fun_print_type = std::function<std::string(const std::vector<std::string>&)>;

int main() {
    expression<gdual_d> ex(1,1,1,6,6,2,function_set<gdual_d>({"sum", "diff", "div", "mul"})(), 67);
    gdual_d x(1.2, "x", 5);
    stream(std::cout, ex({x}), "\n");
}
