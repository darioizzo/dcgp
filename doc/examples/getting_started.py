from dcgpy import expression_gdual_double as expression
from dcgpy import kernel_set_gdual_double as kernel_set
from pyaudi import gdual_double as gdual

# 1- Instantiate a random expression using the 4 basic arithmetic operations
ks = kernel_set(["sum", "diff", "div", "mul"])
ex = expression(1, 1, 1, 6, 6, 2, ks(), 4232123212)

# 2 - Define the symbol set (in our case, 1 input variable named "x") and print the expression
in_sym = ["x"]
print("Expression:", ex(in_sym)[0])

# 3 - Print the simplified expression
print("Simplified expression:", ex.simplify(in_sym))

# 4 - Visualize the dCGP graph
ex.visualize(in_sym)

# 5 - Define a gdual number of value 1.2 and truncation order 2
x = gdual(1.2, "x", 2)

# 6 - Compute the output of the expression and its second derivative in x = 1.2 and print
print("Expression in x=1.2:", ex([x])[0])
print("Second derivative:", ex([x])[0].get_derivative([2]))

# 5 - Mutate the expression with 2 random mutations of active genes and print
ex.mutate_active(2)
print("Mutated expression:", ex(in_sym)[0])
