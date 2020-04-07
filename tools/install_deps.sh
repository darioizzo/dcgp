#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${DCGP_BUILD}" != manylinux* ]]; then
    if [[ "${DCGP_BUILD}" == "OSX*" ]]; then
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
    else
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    fi
    export deps_dir=$HOME/local
    export PATH="$HOME/miniconda/bin:$PATH"
    bash miniconda.sh -b -p $HOME/miniconda
    conda config --add channels conda-forge --force

    conda_pkgs="cmake eigen boost-cpp tbb-devel pagmo-devel audi symengine obake-devel cxx-compiler"

    if [[ "${DCGP_BUILD}" == "Python37" || "${DCGP_BUILD}" == "OSXPython37" ]]; then
        conda_pkgs="$conda_pkgs python=3.7 pyaudi pygmo"
    fi

    # We create the conda environment and activate it
    conda create -q -p $deps_dir -y $conda_pkgs
    source activate $deps_dir

    export deps_dir=$HOME/local
    export PATH="$HOME/miniconda/bin:$PATH"
    export PATH="$deps_dir/bin:$PATH"

    export CXX=clang++
    export CC=clang

    # For python builds, we install pybind11 from the specific commit
    # needed to guarantee interoperability with pyaudi/pygmo
    if [[ "${DCGP_BUILD}" == "Python37" || "${DCGP_BUILD}" == "OSXPython37" ]]; then
        export DCGPY_BUILD_DIR=`pwd`
        git clone https://github.com/pybind/pybind11.git
        cd pybind11
        git checkout 4f72ef846fe8453596230ac285eeaa0ce3278bb4
        mkdir build
        cd build
        cmake \
            -DPYBIND11_TEST=NO \
            -DCMAKE_INSTALL_PREFIX=$DCGPY_BUILD_DIR \
            -DCMAKE_PREFIX_PATH=$DCGPY_BUILD_DIR \
            -DPYTHON_EXECUTABLE=$HOME/local/bin/python3 \
            ..
        make install
        cd ../..
    fi
fi

set +e
set +x
