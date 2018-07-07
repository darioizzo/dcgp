#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

if [[ "${DCGP_BUILD}" == manylinux* ]]; then
    cd ..;
    docker pull ${DOCKER_IMAGE};
    docker run --rm -e TWINE_PASSWORD -e DCGP_BUILD -e TRAVIS_TAG -v `pwd`:/dcgp $DOCKER_IMAGE bash /dcgp/tools/install_docker.sh
fi

set +e
set +x
