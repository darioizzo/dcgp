#include <iostream>
#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>
#include <audi/gdual.hpp>

using namespace dcgp;

int main() {
    // 1- Instantiate a random expression using the 4 basic arithmetic operations
    kernel_set<audi::gdual_d> basic_set({"sum", "diff", "div", "mul"});
    expression<audi::gdual_d> ex(1, 1, 1, 6, 6, 2, basic_set(), 1u);

    // 2 - Define the symbol set (in our case, 1 input variable named "x") and visualize the expression
    std::vector<std::string> in_sym({"x"});
    audi::print("Expression: ", ex(in_sym)[0], "\n");

    // 3 - Define a gdual number of value 1.2 and truncation order 2
    audi::gdual_d x(1.2, "x", 2);

    // 4 - Compute the output of the expression and its second derivative in x = 1.2 and visualize
    audi::print("\nExpression in x=1.2: \n", ex({x})[0], "\n");
    audi::print("\nSecond derivative: ", ex({x})[0].get_derivative({2}), "\n");

    // 5 - Mutate the expression with 2 random mutations of active genes and visualize
    ex.mutate_active(2);
    audi::print("\nMutated expression: ", ex(in_sym)[0], "\n");
}
