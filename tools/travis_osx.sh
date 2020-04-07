#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge --force

conda_pkgs="cmake eigen boost-cpp tbb-devel pagmo-devel audi symengine obake-devel cxx-compiler"

if [[ "${DCGP_BUILD}" == *Python37* ]]; then
    build_cpp_tests=`no`
    conda_pkgs="$conda_pkgs python=3.7 pyaudi pygmo"
else 
    build_cpp_tests=`yes`
fi

conda create -q -p $deps_dir -y $conda_pkgs
source activate $deps_dir

export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
export PATH="$deps_dir/bin:$PATH"

export CXX=clang++
export CC=clang

# Install Pybind11 (making sure its the same used in our pipeline)
export PYAUDI_BUILD_DIR=`pwd`
git clone https://github.com/pybind/pybind11.git
cd pybind11
git checkout 4f72ef846fe8453596230ac285eeaa0ce3278bb4
mkdir build
cd build
pwd
cmake \
    -DPYBIND11_TEST=NO \
    -DCMAKE_INSTALL_PREFIX=$PYAUDI_BUILD_DIR \
    -DCMAKE_PREFIX_PATH=$PYAUDI_BUILD_DIR \
    ..
make install
cd ../..

# Install dcgp
cmake \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_BUILD_TYPE=${DCGP_BUILD_TYPE} \
    -DDCGP_BUILD_DCGP=yes \
    -DDCGP_BUILD_TESTS=$build_cpp_tests \
    -DDCGP_BUILD_EXAMPLES=no \
    ..
make -j2 VERBOSE=1
make install

if [[ "${DCGP_BUILD}" != *Python* ]]; then
    # We run the cpp tests
    ctest -j4 -V
else
    # Else we install and test dcgpy.
    cmake \
        -DCMAKE_INSTALL_PREFIX=$deps_dir \
        -DCMAKE_PREFIX_PATH=$deps_dir \
        -DBoost_NO_BOOST_CMAKE=ON \
        -DCMAKE_BUILD_TYPE=${DCGP_BUILD_TYPE} \
        -DDCGP_BUILD_DCGP=no \
        -DDCGP_BUILD_DCGPY=yes \
        -Dpybind11_DIR=$DCGPY_BUILD_DIR/share/cmake/pybind11 \
        ..
    make -j2 VERBOSE=1
    make install
    # Move out of the build dir.
    cd ../tools
    python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()";
fi

set +e

set +x
