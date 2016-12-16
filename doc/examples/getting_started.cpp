#include <iostream>
#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>

using namespace dcgp;

int main() {
    // 1- Instantiate a random expression using the 4 basic arithmetic operations
    kernel_set<gdual_d> ks({"sum", "diff", "div", "mul"});
    expression<gdual_d> ex(1, 1, 1, 6, 6, 2, ks(), 123212321);

    // 2 - Define the symbol set (in our case, 1 input variable named "x") and visualize the expression
    std::vector<std::string> in_sym({"x"});
    stream(std::cout, "Expression: ", ex(in_sym)[0], "\n");

    // 3 - Define a gdual number of value 1.2 and truncation order 2
    gdual_d x(1.2, "x", 2);

    // 4 - Compute the output of the expression and its second derivative in x = 1.2 and visualize
    stream(std::cout, "Expression in x=1.2: ", ex({x})[0], "\n");
    stream(std::cout, "Second derivative: ", ex({x})[0].get_derivative({2}), "\n");

    // 5 - Mutate the expression with 2 random mutations of active genes and visualize
    ex.mutate_active(2);
    stream(std::cout, "Mutated expression: ", ex(in_sym)[0], "\n");
}
