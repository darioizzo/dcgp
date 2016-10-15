#include <string>

#include "docstrings.hpp"

namespace dcgpy
{

std::string kernel_init_doc(const std::string &type)
{
    return R"(Construct a kernel function from callables.

Args:
    callable_f (``callable: List[)" + type + R"(] -> )" + type + R"(``): a callable taking a list of )" + type + R"( as inputs and returning a )" + type + R"( (the value of the kernel function evaluated on the inputs)
    callable_s (``callable: List[string] -> string``): a callable taking a list of string as inputs and returning a string (the symbolic representation of the kernel function evaluated on the input symbols)
    name (``string``): name of the kernel

Examples:
>>> from dcgpy import *
>>> def my_sum(x):
...     return sum(x)
>>> def print_my_sum(x)
...     s = "+"
...     return "(" + s.join(x) + ") "
>>> my_kernel = kernel_)" + type + R"((my_sum, print_my_sum, "my_sum")
    )";
}

std::string expression_init_doc(const std::string &type)
{
    return R"(A CGP operating on floats
Args:
    inputs (``int``): number of inputs
    outputs (``int``): number of outputs
    rows (``int``): number of rows in the cartesian program
    columns (``int``): number of columns in the cartesian program
    levels_back (``int``): number of levels-back in the cartesian program
    arity (``int``): arity of the kernels
    kernels (``List[pycgp.kernel]``): kernel functions
    seed (``int``): random seed to generate mutations and chromosomes

Examples:
>>> from dcgpy import *
>>> dcgp = expression_)" + type + R"((1,1,1,10,11,2,kernel_set(["sum","diff","mul","div"])(), 32u)
>>> print(dcgp)
...
>>> num_out = dcgp([in])
>>> sym_out = dcgp(["x"])
    )";
}

std::string kernel_set_init_doc(const std::string &type)
{
    return R"(Constructs a set of common kernel functions from their common name. The kernel
functions can be then retrieved via the call operator.

Args:
    kernels (``List[string]``): a list of strings indicating names of kernels to use.
    Available are: "sum", "diff", "mul", "div", "sig", "sin", "log", "exp"

Examples:
>>> from dcgpy import *
>>> kernels = kernel_set_)" + type + R"((["sum", "diff", "mul", "div"])
>>> kernels()[0](["x", "y"])
    )";
}


} // namespace
