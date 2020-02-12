#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

AUDI_VERSION="1.7"
PAGMO_VERSION="2.11.4"

if [[ ${DCGP_BUILD} == *37 ]]; then
	PYTHON_DIR="cp37-cp37m"
	BOOST_PYTHON_LIBRARY_NAME="libboost_python37.so"
	PYTHON_VERSION="37"
	PYTHON_VERSION_DOTTED="3.7"
elif [[ ${DCGP_BUILD} == *36 ]]; then
	PYTHON_DIR="cp36-cp36m"
	BOOST_PYTHON_LIBRARY_NAME="libboost_python36.so"
	PYTHON_VERSION="36"
	PYTHON_VERSION_DOTTED="3.6"
elif [[ ${DCGP_BUILD} == *27mu ]]; then
	PYTHON_DIR="cp27-cp27mu"
	BOOST_PYTHON_LIBRARY_NAME="libboost_python27mu.so"
	PYTHON_VERSION="27"
	PYTHON_VERSION_DOTTED="2.7"
elif [[ ${DCGP_BUILD} == *27 ]]; then
	PYTHON_DIR="cp27-cp27m"
	BOOST_PYTHON_LIBRARY_NAME="libboost_python27.so"
	PYTHON_VERSION="27"
	PYTHON_VERSION_DOTTED="2.7"
else
	echo "Invalid build type: ${DCGP_BUILD}"
	exit 1
fi

# HACK: for python 3.x, the include directory
# is called 'python3.xm' rather than just 'python3.x'.
# This confuses the build system of Boost.Python, thus
# we create a symlink to 'python3.x'.
cd /opt/python/${PYTHON_DIR}/include
PY_INCLUDE_DIR_NAME=`ls`
# If the include dir ends with 'm', create a symlink
# without the 'm'.
if [[ $PY_INCLUDE_DIR_NAME == *m ]]; then
	ln -s $PY_INCLUDE_DIR_NAME `echo $PY_INCLUDE_DIR_NAME|sed 's/.$//'`
fi

cd
cd install

# Python deps
/opt/python/${PYTHON_DIR}/bin/pip install numpy cloudpickle

# Install pybind11
curl -L https://github.com/pybind/pybind11/archive/v2.4.3.tar.gz > v2.4.3
tar xvf v2.4.3 > /dev/null 2>&1
cd pybind11-2.4.3
mkdir build
cd build
cmake ../ -DPYBIND11_TEST=OFF > /dev/null
make install > /dev/null 2>&1
cd ../..

# Install pagmo and pygmo
curl -L  https://github.com/esa/pagmo2/archive/v${PAGMO_VERSION}.tar.gz > pagmo2.tar.gz
tar xzf pagmo2.tar.gz
cd pagmo2-${PAGMO_VERSION}
mkdir build_pagmo
cd build_pagmo
cmake -DBoost_NO_BOOST_CMAKE=ON \
	-DPAGMO_WITH_EIGEN3=yes \
	-DPAGMO_WITH_NLOPT=yes \
	-DPAGMO_WITH_IPOPT=yes \
	-DCMAKE_BUILD_TYPE=Release ../;
make -j2 install
cd ../
mkdir build_pygmo
cd build_pygmo
cmake -DBoost_NO_BOOST_CMAKE=ON \
	-DCMAKE_BUILD_TYPE=Release \
	-DPAGMO_BUILD_PYGMO=yes \
	-DPAGMO_BUILD_PAGMO=no \
	-DBoost_PYTHON${PYTHON_VERSION}_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIBRARY_NAME} \
	-DPYTHON_EXECUTABLE=/opt/python/${PYTHON_DIR}/bin/python \
	-DYACMA_PYTHON_MODULES_INSTALL_PATH=/opt/python/${PYTHON_DIR}/lib/python${PYTHON_VERSION_DOTTED}/site-packages ../;
make -j2 install
cd ../..

# Install audi
curl -L https://github.com/darioizzo/audi/archive/v${AUDI_VERSION}.tar.gz > v${AUDI_VERSION}
tar xvf v${AUDI_VERSION} > /dev/null 2>&1
cd audi-${AUDI_VERSION}
mkdir build
cd build
cmake -DBoost_NO_BOOST_CMAKE=ON -DAUDI_BUILD_AUDI=yes -DAUDI_BUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release ../
make install > /dev/null 2>&1
cd ../..

# Install dcgp headers
cd /dcgp
mkdir build_dcgp
cd build_dcgp
cmake -DBoost_NO_BOOST_CMAKE=ON -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release ../
make install

# Compile and install dcgpy (build directory is created by .travis.yml)
cd /dcgp
cd build
cmake -DBoost_NO_BOOST_CMAKE=ON \
	-DCMAKE_BUILD_TYPE=Release \
	-DDCGP_BUILD_DCGP=no \
	-DDCGP_BUILD_DCGPY=yes \
	-DBoost_PYTHON${PYTHON_VERSION}_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIBRARY_NAME} \
	-DPYTHON_EXECUTABLE=/opt/python/${PYTHON_DIR}/bin/python ../;
make -j2 install


# Compile wheels
cd wheel
# Copy the installed pyaudi files, wherever they might be in /usr/local,
# into the current dir.
cp -a `find /usr/local/lib -type d -iname 'dcgpy'` ./
# Create the wheel and repair it.
/opt/python/${PYTHON_DIR}/bin/python setup.py bdist_wheel
auditwheel repair dist/dcgpy* -w ./dist2
# Try to install it and run the tests.
cd /
/opt/python/${PYTHON_DIR}/bin/pip install /dcgp/build/wheel/dist2/dcgpy*
/opt/python/${PYTHON_DIR}/bin/python -c "from dcgpy import test; test.run_test_suite(); import pygmo; pygmo.mp_island.shutdown_pool(); pygmo.mp_bfe.shutdown_pool()"

# Upload in PyPi
# This variable will contain something if this is a tagged build (vx.y.z), otherwise it will be empty.
export DCGP_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`
if [[ "${DCGP_RELEASE_VERSION}" != "" ]]; then
    echo "Release build detected, uploading to PyPi."
    /opt/python/${PYTHON_DIR}/bin/pip install twine
    /opt/python/${PYTHON_DIR}/bin/twine upload -u darioizzo /dcgp/build/wheel/dist2/dcgpy*
fi
