#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Install conda+deps.
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda_pkgs="cmake boost-cpp eigen pagmo-devel tbb audi symengine obake-devel cxx-compiler pybind11 pybind11-abi python=3.8 sympy matplotlib pygmo pyaudi"
conda create -y -q -p $deps_dir $conda_pkgs
source activate $deps_dir

# Create the build dir and cd into it.
mkdir build
cd build

# Clang build
CXX=clang++ CC=clang cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DDCGP_BUILD_TESTS=no \
    -DDCGP_BUILD_EXAMPLES=no \
    ..
make install

cd ..
mkdir build_python
cd build_python

# Now the python bindings.
CXX=clang++ CC=clang cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DDCGP_BUILD_DCGP=no \
    -DDCGP_BUILD_DCGPY=yes \
    ..

make -j2 VERBOSE=1
make install

python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"


set +e
set +x