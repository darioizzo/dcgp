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

int main() {
    basis_function<double> fun1(my_sum<double>, print_my_sum, "sum");
    basis_function<gdual_d> fun2(my_sum<gdual_d>, print_my_sum, "sum");
    basis_function<gdual_v> fun3(my_sum<gdual_v>, print_my_sum, "sum");
    std::cout << fun1({1.,2.,3.}) << "\n";
    std::cout << fun2({gdual_d(1.),gdual_d(2.),gdual_d(3.)}) << "\n";
    std::cout << fun3({gdual_v(std::vector<double>{1.,2.}), gdual_v(std::vector<double>{2.,3.}), gdual_v(std::vector<double>{3.,4.})}) << "\n";
    function_set<double> basic_set1({"sum","diff","mul"});
    expression<double> my_expr1(1u, 1u, 1u, 6u, 6u, 2u, basic_set1(), 23u);
    stream(std::cout, my_expr1(std::vector<std::string>{"x"}), "\n");
    stream(std::cout, my_expr1(std::vector<double>{1.}), "\n");
    function_set<gdual_d> basic_set2({"sum","diff","mul"});
    expression<gdual_d> my_expr2(1u, 1u, 1u, 6u, 6u, 2u, basic_set2(), 23u);
    stream(std::cout, my_expr2(std::vector<std::string>{"x"}), "\n");
    stream(std::cout, my_expr2(std::vector<gdual_d>{gdual_d(1.)}), "\n");
    function_set<gdual_v> basic_set3({"sum","diff","mul"});
    expression<gdual_v> my_expr3(1u, 1u, 1u, 6u, 6u, 2u, basic_set3(), 23u);
    stream(std::cout, my_expr3(std::vector<std::string>{"x"}), "\n");
    stream(std::cout, my_expr3(std::vector<gdual_v>{gdual_v(1.)}), "\n");
}
