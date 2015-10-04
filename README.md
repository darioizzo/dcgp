# d-CGP
Implementation of the differential Cartesian Genetic Programming (d-CGP)

The d-CGP is a recent development in the field of Genetic Programming that adds the information about the derivatives of the output nodes (the programs, or expressions encoded) with respect to the input nodes (the input values) and weights. In doing so, it enables a number of new applications currently the subject of active research.

 * The evolution of the genetic program can now be helped by using the information on the derivatives, enabling for the equivalent of backpropagation in Neural Networks.
 * The fitness function can be defined in terms of the derivatives, allowing to go beyond simple regression tasks and, instead, solve differential equations, learn differential models, capture conserved quantities in dynamical systems.

### Dependencies
Several and tested dependencies are necessary to succesfully use d-CGP
 * Audi, headers only library - (git clone https://github.com/darioizzo/audi.git)
 * Piranha, headers only library - (git clone https://github.com/bluescarni/piranha.git)
 * Boost, headers only

### Comparison to the CGP-Library
If all below statements are true:
 * You do not care about knowing derivatives of your encoded function
 * You do not care about run-time capabilities
 * You do not care about the Python module
 * You do not care about the possibility of defining your kernel functions as complex functors (e.g. CGP expressions.)
 * You do not care about thread-safety

then you should consider using, instead, Andrew Turner's CGP-Library (http://www.cgplibrary.co.uk/files2/About-txt.html) which is, roughly, twice as fast to compute a CGP expression as it makes use of function pointers rather than a std::function to define the kernel functions.
