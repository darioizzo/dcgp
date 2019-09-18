#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${DCGP_BUILD}" != manylinux* ]]; then
    export deps_dir=$HOME/local
    export PATH="$HOME/miniconda/bin:$PATH"
    export PATH="$deps_dir/bin:$PATH"
fi

if [[ "${DCGP_BUILD}" == "ReleaseGCC" ]]; then
    cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Release -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes -DDCGP_BUILD_EXAMPLES=no -DCMAKE_CXX_FLAGS="-fuse-ld=gold" ../;
    make -j2 VERBOSE=1;
    ctest -VV;
elif [[ "${DCGP_BUILD}" == "DebugGCC" ]]; then
    cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Debug -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes -DDCGP_BUILD_EXAMPLES=no -DCMAKE_CXX_FLAGS="-fsanitize=address -fuse-ld=gold" ../;
    make -j2 VERBOSE=1;
    ctest -VV;
elif [[ "${DCGP_BUILD}" == "CoverageGCC" ]]; then
    cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Debug -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes -DDCGP_BUILD_EXAMPLES=no -DCMAKE_CXX_FLAGS="--coverage -fuse-ld=gold" ../;
    make -j2 VERBOSE=1;
    ctest -VV;
    bash <(curl -s https://codecov.io/bash) -x gcov-5;
elif [[ "${PAGMO_PLUGINS_NONFREE_BUILD}" == "DebugClang" ]]; then
    CXX=clang++-7 CC=clang-7 cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Debug --DCMAKE_BUILD_TYPE=Debug -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes ../;
    make -j2 VERBOSE=1;
    ctest -VV;
elif [[ "${DCGP_BUILD}" == "ReleaseClang" ]]; then
    CXX=clang++-7 CC=clang-7 cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_TYPE=Debug -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes ../;
    make -j2 VERBOSE=1;
    ctest -VV;
elif [[ "${DCGP_BUILD}" == "OSXDebug" ]]; then
    CXX=clang++ CC=clang cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_BUILD_TYPE=Debug -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes -DDCGP_BUILD_EXAMPLES=no -DCMAKE_CXX_FLAGS="-g0 -O2" ../;
    make -j2 VERBOSE=1;
    ctest -VV;
elif [[ "${DCGP_BUILD}" == "OSXRelease" ]]; then
    CXX=clang++ CC=clang cmake -DCMAKE_PREFIX_PATH=$deps_dir -DBoost_NO_BOOST_CMAKE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_TYPE=Debug -DDCGP_BUILD_DCGP=yes -DDCGP_BUILD_TESTS=yes -DDCGP_BUILD_EXAMPLES=no ../;
    make -j2 VERBOSE=1;
    ctest -VV;
elif [[ "${DCGP_BUILD}" == manylinux* ]]; then
    cd ..;
    docker pull ${DOCKER_IMAGE};
    docker run --rm -e TWINE_PASSWORD -e DCGP_BUILD -e TRAVIS_TAG -v `pwd`:/dcgp $DOCKER_IMAGE bash /dcgp/tools/install_docker.sh
fi

set +e
set +x
