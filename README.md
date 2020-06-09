[![Travis (.org)](https://img.shields.io/travis/darioizzo/dcgp?style=for-the-badge&logo=travis)](https://travis-ci.org/darioizzo/dcgp)
[![Azure DevOps builds](https://img.shields.io/azure-devops/build/darioizzo/c3218f2b-0cd5-476a-af68-5abe1e6a7c14/2?style=for-the-badge)](https://dev.azure.com/darioizzo/dcgp/_build)
[![PyPI](https://img.shields.io/pypi/v/dcgpy?style=for-the-badge)](https://pypi.python.org/pypi?:action=display&name=dcgpy&version=1.0.1)
[![Conda (channel only)](https://img.shields.io/conda/vn/conda-forge/dcgp-python?style=for-the-badge)](https://anaconda.org/conda-forge/dcgp-python)

[![Gitter](https://img.shields.io/gitter/room/esa/dcgp?logo=gitter-white&style=for-the-badge)](https://gitter.im/esa/dcgp)
![GitHub commit activity](https://img.shields.io/github/commit-activity/y/darioizzo/dcgp?style=for-the-badge)

[![DOI](https://zenodo.org/badge/38923383.svg)](https://zenodo.org/badge/latestdoi/38923383)

# dCGP
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

Installation guide
==================

C++
---

dCGP is a header-only library which has the following third party dependencies:

* [Boost](http://www.boost.org/), various C++ utilities (>=1.72).
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), linear algebra library (>=3.3.0)
* [Pagmo](https://github.com/esa/pagmo2), parallel optimization library (>=2.15).
* [tbb](https://github.com/intel/tbb), lets you easily write parallel C++ programs that take full advantage of multicore performance (>=2020.1).
* [AuDi](http://darioizzo.github.io/audi/), high order automated differentiation library (>=1.8).
* [Symengine](https://github.com/symengine/symengine), symbolic manipulation of math expressions (>=0.6).

After making sure the dependencies above are installed and found in your system, you may download
the latest dCGP code via git:

```bash
git clone https://github.com/darioizzo/dcgp.git
```
and configure your build using CMake. 

When done, type (in your build directory):

```bash
 make 
```
When finished, to run the tests type:

```bash
make test
```

If succesfull, you may now install cgp:

```bash
make install
```

The headers will be installed in the CMAKE_INSTALL_PREFIX/include directory. 

Python
------
The main functionalities of dCGP are exposed into a Python module called ``dcgpy`` which
can be installed from conda (OSx, linux and Win), pip (only linux) or by building the module.

### Installing with conda

``dcgpy`` is available in the `conda <https://conda.io/en/latest/>`__ package manager
from the `conda-forge <https://conda-forge.org/>`__ channel. A single package is available:

* `dcgp-python <https://anaconda.org/conda-forge/dcccgp-python>`__, which contains the ``dcgpy`` python module.

In order to install ``dcgpy`` via conda, you just need
to add ``conda-forge`` to the channels:

```
conda config --add channels conda-forge
conda install dcgp-python
```

Please refer to the `conda documentation <https://conda.io/en/latest/>`__ for instructions
on how to setup and manage your conda installation.

You may test the successfull installation by running the python tests typing:

```
python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"
```

### Installing with pip (deprecated)

We also provide the pip packages (mainly for linux 64 bit architectures and older versions).
Check on the `PyPi dcgpy page <https://pypi.org/project/dcgpy/>`_ if the needed package is provided.

```
pip install dcgpy
```

### Building


To build the module you need to have the Boost Python libraries installed and to activate the BUILD_DCGPY option from within CMake (and deselect BUILD_DCGP)

Check carefully what Python version is detected and what libraries are linked to. In particular select the correct boost_python
according to the Python version (2 or 3) you want to compile the module for.

The CMAKE_INSTALL_PREFIX will be used to construct the final location of headers and Python module after install.

When done, type (in your build directory):

```
make install
```

To check that all went well fire-up your Python console and try the example in :ref:`quick-start example <getting_started_py>`.


## Other CGP libraries
### Comparison to the CGP-Library
If all below statements are true:
 * You do not care about knowing derivatives of your encoded program
 * You do not care about run-time capabilities
 * You do not care about the Python module
 * You do not care about the possibility of defining your kernel functions as complex functors (e.g. CGP expressions.)
 * You do not care about thread-safety

then you should consider using, instead, Andrew Turner's CGP-Library (http://www.cgplibrary.co.uk/files2/About-txt.html) which is, roughly, twice as fast to compute a CGP expression as it makes use of function pointers rather than a type-erased function wrapper to define the kernel functions.
