#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${DCGP_BUILD}" != manylinux* ]]; then
    if [[ "${TRAVIS_OS_NAME}" == "osx" ]]; then
        wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;
    else
        wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
    fi
    export deps_dir=$HOME/local
    export PATH="$HOME/miniconda/bin:$PATH"
    bash miniconda.sh -b -p $HOME/miniconda
    conda config --add channels conda-forge --force

    # obake-devel is needed as far as the conda package audi does not list it as a dependency
    conda_pkgs="cmake eigen boost-cpp tbb-devel pagmo-devel audi symengine obake-devel cxx-compiler"

    if [[ "${DCGP_BUILD}" == "Python37" || "${DCGP_BUILD}" == "OSXPython37" ]]; then
        conda_pkgs="$conda_pkgs python=3.7 pyaudi pygmo"
    fi

    # We create the conda environment and activate it
    conda create -q -p $deps_dir -y
    source activate $deps_dir
    conda install $conda_pkgs -y

    # We install pybind11 from the specific commit needed to guarantee interoperability with pyaudi/pygmo
    export DCGPY_BUILD_DIR=`pwd`
    git clone https://github.com/pybind/pybind11.git
    cd pybind11
    git checkout 4f72ef846fe8453596230ac285eeaa0ce3278bb4
    mkdir build
    cd build

    python -c "from distutils import sysconfig as s;import sys;import struct;
        print('.'.join(str(v) for v in sys.version_info));
        print(sys.prefix);
        print(s.get_python_inc(plat_specific=True));
        print(s.get_python_lib(plat_specific=True));
        print(s.get_config_var('SO'));
        print(hasattr(sys, 'gettotalrefcount')+0);
        print(struct.calcsize('@P'));
        print(s.get_config_var('LDVERSION') or s.get_config_var('VERSION'));
        print(s.get_config_var('LIBDIR') or '');
        print(s.get_config_var('MULTIARCH') or '');
        "

    cmake \
        -DPYBIND11_TEST=NO \
        -DCMAKE_INSTALL_PREFIX=$DCGPY_BUILD_DIR \
        -DCMAKE_PREFIX_PATH=$DCGPY_BUILD_DIR \
        -DPYTHON_EXECUTABLE=$HOME/miniconda/bin/python3 \
        -DPYTHON_LIBRARY=$HOME/miniconda/lib/libpython3.so \
        -DPYBIND11_PYTHON_VERSION=
        ..
    make install
    cd ../..
fi

set +e
set +x
