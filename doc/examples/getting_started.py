from pyaudi import gdual_double as gdual
from pyaudi import exp, log, cbrt

# We want to compute the Taylor expansion of a function f (and thus all derivatives) at x=2, y=3
# 1 - Define the generalized dual numbers (7 is the truncation order, i.e. the maximum order of derivation we will need)

x = gdual(2, "x", 7);
y = gdual(3, "y", 7);

# 2 - Compute your function as usual
f = exp(x*x + cbrt(y) / log(x*y));

# 3 - Inspect the results (this does not require any more computations)
print("Taylor polynomial: " + str(f))                          # This is the Taylor expansion of f (truncated at the 7th order)
print("Derivative value [1,0]: " + str(f.get_derivative([1,0])))     # This is the value of the derivative (d / dx)
print("Derivative value [4,3]: " + str(f.get_derivative([4,3])))     # This is the value of the mixed derivative (d^7 / dx^4dy^3)

# 4 - Using the dictionary interface (note the presence of the "d" before all variables)
print("Derivative value [1,0]: " + str(f.get_derivative({"dx":1})))     # This is the value of the derivative (d / dx)
print("Derivative value [4,3]: " + str(f.get_derivative({"dx":4, "dy":3})))     # This is the value of the mixed derivative (d^7 / dx^4dy^3)
