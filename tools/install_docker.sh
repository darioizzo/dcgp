#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

CMAKE_VERSION="3.11.1"
EIGEN3_VERSION="3.3.4"
BOOST_VERSION="1.67.0"
AUDI_VERSION="1.5"
PIRANHA_VERSION="0.11"


if [[ ${DCGP_BUILD} == *37 ]]; then
	PYTHON_DIR="cp37-cp37m"
elif [[ ${DCGP_BUILD} == *36 ]]; then
	PYTHON_DIR="cp36-cp36m"
elif [[ ${DCGP_BUILD} == *35 ]]; then
	PYTHON_DIR="cp35-cp35m"
elif [[ ${DCGP_BUILD} == *27 ]]; then
	PYTHON_DIR="cp27-cp27mu"
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
mkdir install
cd install

# Install Boost
curl -L http://dl.bintray.com/boostorg/release/${BOOST_VERSION}/source/boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2 > boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2
tar xjf boost_`echo ${BOOST_VERSION}|tr "." "_"`.tar.bz2
cd boost_`echo ${BOOST_VERSION}|tr "." "_"`
sh bootstrap.sh --with-python=/opt/python/${PYTHON_DIR}/bin/python > /dev/null
./bjam --toolset=gcc link=shared threading=multi cxxflags="-std=c++11" variant=release --with-python --with-serialization --with-iostreams --with-regex --with-chrono --with-timer --with-test --with-system -j2 install > /dev/null
cd ..

# Install gmp (before mpfr as its used by it)
curl -L https://gmplib.org/download/gmp/gmp-6.1.2.tar.bz2 > gmp-6.1.2.tar.bz2
tar xvf gmp-6.1.2.tar.bz2  > /dev/null 2>&1
cd gmp-6.1.2 > /dev/null 2>&1
./configure --enable-fat > /dev/null 2>&1
make > /dev/null 2>&1
make install > /dev/null 2>&1
cd ..

# Install mpfr
curl -L http://www.mpfr.org/mpfr-3.1.6/mpfr-3.1.6.tar.gz > mpfr-3.1.6.tar.gz
tar xvf mpfr-3.1.6.tar.gz > /dev/null 2>&1
cd mpfr-3.1.6
./configure > /dev/null 2>&1
make > /dev/null 2>&1
make install > /dev/null 2>&1
cd ..

# Install CMake
curl -L https://github.com/Kitware/CMake/archive/v${CMAKE_VERSION}.tar.gz > v${CMAKE_VERSION}
tar xzf v${CMAKE_VERSION} > /dev/null 2>&1
cd CMake-${CMAKE_VERSION}/
./configure > /dev/null
gmake -j2 > /dev/null
gmake install > /dev/null
cd ..

# Install Eigen
curl -L http://bitbucket.org/eigen/eigen/get/${EIGEN3_VERSION}.tar.gz > ${EIGEN3_VERSION}
tar xzf ${EIGEN3_VERSION} > /dev/null 2>&1
cd eigen*
mkdir build
cd build
cmake ../ > /dev/null
make install > /dev/null
cd ..
cd ..

# Install piranha
curl -L https://github.com/bluescarni/piranha/archive/v${PIRANHA_VERSION}.tar.gz > v${PIRANHA_VERSION}
tar xvf v${PIRANHA_VERSION} > /dev/null 2>&1
cd piranha-${PIRANHA_VERSION}
mkdir build
cd build
cmake ../ > /dev/null
make install > /dev/null 2>&1
cd ..

# Install audi
git clone https://github.com/darioizzo/audi.git > /dev/null 2>&1
cd audi
mkdir build
cd build
cmake -DBUILD_TESTS=no ../
make install > /dev/null 2>&1

# Install and compile dcgpy
cd /dcgp
mkdir build
cd build
# The include directory for py3 is X.Xm, while for py2 is X.X
if [[ "${PYTHON_VERSION}" != "2.7" ]]; then
    cmake -DBUILD_DCGPY=yes -DBUILD_TESTS=no -DBUILD_MAIN=no -DCMAKE_INSTALL_PREFIX=/dcgp/local -DCMAKE_BUILD_TYPE=Release -DBoost_PYTHON_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIB_NAME} -DPYTHON_INCLUDE_DIR=${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}m/ -DPYTHON_EXECUTABLE=${PATH_TO_PYTHON}/bin/python  ../
else
    cmake -DBUILD_DCGPY=yes -DBUILD_TESTS=no -DBUILD_MAIN=no -DCMAKE_INSTALL_PREFIX=/dcgp/local -DCMAKE_BUILD_TYPE=Release -DBoost_PYTHON_LIBRARY_RELEASE=/usr/local/lib/${BOOST_PYTHON_LIB_NAME} -DPYTHON_INCLUDE_DIR=${PATH_TO_PYTHON}/include/python${PYTHON_VERSION}/ -DPYTHON_EXECUTABLE=${PATH_TO_PYTHON}/bin/python  ../
fi
make
make install

# Compile wheels
cd /dcgp/build/wheel
cp -R /dcgp/local/lib/python${PYTHON_VERSION}/site-packages/dcgpy ./
# The following line is needed as a workaround to the auditwheel problem KeyError = .lib
# Using and compiling a null extension module (see manylinux_wheel_setup.py)
# fixes the issue (TODO: probably better ways?)
touch dummy.cpp

# We install required dependencies (do it here, do not let pip install do it)
${PATH_TO_PYTHON}/bin/pip install numpy pyaudi
${PATH_TO_PYTHON}/bin/pip wheel ./ -w wheelhouse/
# Bundle external shared libraries into the wheels (only py35 has auditwheel)
/opt/python/cp35-cp35m/bin/auditwheel repair wheelhouse/dcgpy*.whl -w ./wheelhouse2/
# Install packages (not sure what --no-index -f does, should also work without, but just in case)
${PATH_TO_PYTHON}/bin/pip install dcgpy --no-index -f wheelhouse2
# Test
cd /
${PATH_TO_PYTHON}/bin/python -c "from dcgpy import test; test.run_test_suite()"

# Upload in PyPi
# This variable will contain something if this is a tagged build (vx.y.z), otherwise it will be empty.
#export DCGP_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`
if [[ "${DCGP_RELEASE_VERSION}" != "" ]]; then
    echo "Release build detected, uploading to PyPi."
    cd dcgp/build/wheel
    ${PATH_TO_PYTHON}/bin/pip install twine
    ${PATH_TO_PYTHON}/bin/twine upload -u darioizzo wheelhouse2/dcgpy*.whl
fi
