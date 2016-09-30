#include <string>

#include "docstrings.hpp"

namespace dcgpy
{

std::string basis_function_init_doc()
{
    return R"(Construct a kernel function from callables.

Args:
    callable_f (``callable: List[float] -> float``): a callable taking a list of floats as inputs and returning a float (the value of the kernel function evaluated on the inputs)
    callable_s (``callable: List[string] -> string``): a callable taking a list of string as inputs and returning a string (the symbolic representation of the kernel function evaluated on the input symbols)
    name (``string``): name of the kernel

Examples:
>>> from dcgpy import *
>>> def my_sum(x):
...     return sum(x)
>>> def print_my_sum(x)
...     s = "+"
...     return "(" + s.join(x) + ") "
>>> my_kernel = kernel(my_sum, print_my_sum, "my_sum")
    )";
}

std::string expression_init_doc()
{
    return R"(A CGP operating on floats
Args:
    in (``int``): number of inputs
    out (``int``): number of outputs
    rows (``int``): number of rows in the cartesian program
    columns (``int``): number of columns in the cartesian program
    levels-back (``int``): number of levels-back in the cartesian program
    arity (``int``): arity of the kernels
    kernels (``List[pycgp.kernel]``): kernel functions
    seed (``int``): random seed to generate mutations and chromosomes

Examples:
>>> from dcgpy import *
>>> cgp = expression(1,1,1,10,11,2,function_set(["sum","diff","mul","div"])(), 32u)
>>> print(cgp)
...
>>> num_out = cgp([2.1])
>>> sym_out = cgp(["x"])
    )";
}

std::string function_set_init_doc()
{
    return R"(Constructs a set of common kernel functions from their common name. The kernel
functions can be then retrieved via the call operator.

Args:
    kernels (``List[string]``): a list of strings indicating names of kernels to use.
    Available are: "sum", "diff", "mul", "div", "sig", "sin", "log", "exp"

Examples:
>>> from dcgpy import *
>>> kernels = function_set(["sum", "diff", "mul", "div"])
>>> kernels()[0](["x", "y"])
    )";
}


} // namespace
