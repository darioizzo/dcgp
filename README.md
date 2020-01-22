[![Travis (.org)](https://img.shields.io/travis/darioizzo/dcgp?style=for-the-badge&logo=travis)](https://travis-ci.org/darioizzo/dcgp)
[![Azure DevOps builds](https://img.shields.io/azure-devops/build/darioizzo/c3218f2b-0cd5-476a-af68-5abe1e6a7c14/2?style=for-the-badge)](https://dev.azure.com/darioizzo/dcgp/_build)
[![PyPI](https://img.shields.io/pypi/v/dcgpy?style=for-the-badge)](https://pypi.python.org/pypi?:action=display&name=dcgpy&version=1.0.1)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/dcgp-python?style=for-the-badge)](https://anaconda.org/conda-forge/dcgp-python)

[![DOI](https://zenodo.org/badge/38923383.svg)](https://zenodo.org/badge/latestdoi/38923383)

# d-CGP
Implementation of differentiable Cartesian Genetic Programming (dCGP)

The dCGP is a development in the field of Genetic Programming that adds the information about the derivatives of the program output with respect to the input values and various parameters (weights, biases, etc..). In doing so, it enables a number of new applications currently the subject of active research.

 * The representation of deep neural networks using a dCGPANN allows the whole network to be encoded and evolved, including the connection topology, the weights, the biases, etc, .... as well as to learn the weights and biases using backpropagation.
 * The solution to boundary values problems, differential equations etc. can be encoded in a dCGP and evolved against different boundary conditions or initial conditions
 * Prime integrals of motion can be represented by a dCGP and learned
 * Symbolic regression tasks can learn ephemeral constants using backprop or even higher order methods.
 
The first research paper describing dCGP use to solve symbolic regressions problems such is:

Izzo, Dario, Francesco Biscani, and Alessio Mereta. "Differentiable Genetic Programming." arXiv preprint arXiv:1611.04766 (2016).

Preliminary documentation can be found at http://darioizzo.github.io/dcgp/

A web based version of dCGP can be found here: https://esa.github.io/dcgp-web/ thanks to Mike Heddes!

## dcgpy
If you have a win 64bit system or a linux based system (32 or 64 bits), the python package dcgpy (python binding of the C++ code) can be installed via:

 ```pip install dcgpy```

otherwise you will have to compile it by activating the BUILD_DCGPY option in CMake

## Compiling the source code or using the header only library
### Dependencies
Several dependencies are necessary to successfully compile d-CGP
 * Audi, headers only library - (git clone https://github.com/darioizzo/audi.git)
 * Boost, headers only
 
After making sure the dependencies above are installed in your system, you may download the latest dCGP version via git:

```git clone https://github.com/darioizzo/d-CGP.git```

and configure your build using CMake. 

```cd d-CGP```

```mkdir build```

```cd build```

```ccmake ../```

To build the python package activate the BUILD_DCGPY option.

When done, type (in your build directory):

```make install```

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. To check that all went well compile the quick-start example.

## Other CGP libraries
### Comparison to the CGP-Library
If all below statements are true:
 * You do not care about knowing derivatives of your encoded program
 * You do not care about run-time capabilities
 * You do not care about the Python module
 * You do not care about the possibility of defining your kernel functions as complex functors (e.g. CGP expressions.)
 * You do not care about thread-safety

then you should consider using, instead, Andrew Turner's CGP-Library (http://www.cgplibrary.co.uk/files2/About-txt.html) which is, roughly, twice as fast to compute a CGP expression as it makes use of function pointers rather than a type-erased function wrapper to define the kernel functions.
