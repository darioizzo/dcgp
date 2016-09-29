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
    fun_type my_cpp_callable1 = [](std::vector<double> x) { return x[0];};
    fun_print_type my_cpp_callable2 = [](std::vector<std::string> x) { return x[0];};
    auto a = basis_function<double>(my_cpp_callable1, my_cpp_callable2, "giggio");
}
