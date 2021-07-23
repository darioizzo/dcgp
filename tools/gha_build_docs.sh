#!/usr/bin/env bash

# Echo each command
set -x

# Exit on error.
set -e

export deps_dir=$HOME/local
export PATH="$HOME/miniconda/bin:$PATH"
conda_extra_pkgs="sphinx breathe sphinxcontrib-napoleon nbsphinx doxygen sphinx_rtd_theme"
source activate $deps_dir
conda install $conda_extra_pkgs -y

# Move to docs folder
cd doc
# run doxygen
cd doxygen
doxygen
# run sphinx
cd ../sphinx
make html

set +e
set +x