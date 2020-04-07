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

if [[ "${DCGP_BUILD}" == "Python37" || "${DCGP_BUILD}" == "OSXPython37" ]]; then
    conda_pkgs="$conda_pkgs python=3.7 pyaudi pygmo"
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
    -DDCGP_BUILD_TESTS=yes \
    -DDCGP_BUILD_EXAMPLES=no \
    ..
make -j2 VERBOSE=1
make install
ctest -j4 -V

set +e
set +x

