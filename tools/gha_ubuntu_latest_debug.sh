#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

# Core deps.
sudo apt-get install build-essential wget

# Install conda+deps.
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
bash miniconda.sh -b -p $HOME/miniconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda_pkgs="cmake boost-cpp eigen pagmo-devel tbb audi symengine obake-devel cxx-compiler pygmo"
conda create -q -p $deps_dir -y
source activate $deps_dir
conda install $conda_pkgs -y

# Create the build dir and cd into it.
mkdir build
cd build

# GCC build. The -lrt flag seems a necessary addition when using conda compilers installed via cxx-compiler
LDFLAGS="-lrt ${LDFLAGS}"; cmake \
    -DBoost_NO_BOOST_CMAKE=ON \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_INSTALL_PREFIX=$deps_dir \
    -DCMAKE_PREFIX_PATH=$deps_dir \
    -DDCGP_BUILD_TESTS=yes \
    -DDCGP_BUILD_EXAMPLES=no \
    ..
make -j1 VERBOSE=1

python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"

set +e
set +x