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

# Install tbb (headers will be in include, libraries in build/release build/debug)
curl -L https://github.com/01org/tbb/archive/2019.tar.gz >tbb-2019.tar.gz
tar xvf tbb-2019.tar.gz > /dev/null 2>&1
cd tbb-2019
make > /dev/null 2>&1
cd build
mv *_release release
mv *_debug debug
cd ../../

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
curl -L https://github.com/darioizzo/audi/archive/v${AUDI_VERSION}.tar.gz > v${AUDI_VERSION}
tar xvf v${AUDI_VERSION} > /dev/null 2>&1
cd audi-${AUDI_VERSION}
mkdir build
cd build
cmake -DAUDI_BUILD_AUDI=yes -DAUDI_BUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release ../
make install > /dev/null 2>&1

# Python deps
/opt/python/${PYTHON_DIR}/bin/pip install numpy
sleep 20

# Install dcgp headers
cd /dcgp
mkdir build_dcgp
cd build_dcgp
cmake -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=no -DCMAKE_BUILD_TYPE=Release -DTBB_ROOT_DIR=/root/install/tbb-2019 -DTBB_LIBRARY=/root/install/tbb-2019/build/release/ ../
make install

# Compile and install dcgpy (build directory is created by .travis.yml)
cd /dcgp
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DDCGP_BUILD_DCGP=no -DDCGP_BUILD_DCGPY=yes -DPYTHON_EXECUTABLE=/opt/python/${PYTHON_DIR}/bin/python -DTBB_ROOT_DIR=/root/install/tbb-2019 -DTBB_LIBRARY=/root/install/tbb-2019/build/release/ ../;
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
/opt/python/${PYTHON_DIR}/bin/python -c "from dcgpy import test; test.run_test_suite()"

# Upload in PyPi
# This variable will contain something if this is a tagged build (vx.y.z), otherwise it will be empty.
export DCGP_RELEASE_VERSION=`echo "${TRAVIS_TAG}"|grep -E 'v[0-9]+\.[0-9]+.*'|cut -c 2-`
if [[ "${DCGP_RELEASE_VERSION}" != "" ]]; then
    echo "Release build detected, uploading to PyPi."
    /opt/python/${PYTHON_DIR}/bin/pip install twine
    /opt/python/${PYTHON_DIR}/bin/twine upload -u darioizzo /dcgp/build/wheel/dist2/dcgpy*
fi
