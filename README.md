# d-CGP
Implementation of the differential CGP (Cartesian Genetic Programming)

The d-CGP adds to the classical Genetic Programming apporach the information about the derivatives of the output nodes (the programs, or expressions encoded) with respect to the input nodes (the input values and weights). In doing so, it enables a number of new applications currently the subject of active research.

 * The evolution of the genetic program can now be helped by using the information on the derivatives, enabling for the equivalent of backpropagation in Neural Networks
 * The fitness function can be defined in terms of the derivatives, allowing to go beyond simple regression tasks and, insted, solve differential equations, learn differential models, capture conserved quantities in dynamical systems.

### Dependencies
Several 'unusual' and tested dependencies are necessary to succesfully use d-CGP
 * Audi, headers only library - (git clone https://github.com/darioizzo/audi.git)
 * Piranha, headers only library - (git clone https://github.com/bluescarni/piranha.git)
 * Boost, headers only
