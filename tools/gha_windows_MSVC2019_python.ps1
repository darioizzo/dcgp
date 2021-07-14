# Powershell script

# Install conda environment
conda config --set always_yes yes
conda create --name dcgpy cmake boost-cpp eigen pagmo-devel tbb audi symengine obake-devel cxx-compiler pybind11 python=3.8 sympy matplotlib clang ninja
conda activate dcgpy

#build and run the dcgp tests
mkdir build
cd build
cmake `
    -G "Ninja" `
    -DCMAKE_C_COMPILER=clang-cl `
    -DCMAKE_CXX_COMPILER=clang-cl `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\dcgpy `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\dcgpy `
    -DDCGP_BUILD_TESTS=no `
    -DDCGP_BUILD_EXAMPLES=no ..

cmake --build . --target install --config Release

cd ..
mkdir build_python
cd build_python

cmake `
    -G "Ninja" `
    -DCMAKE_C_COMPILER=clang-cl `
    -DCMAKE_CXX_COMPILER=clang-cl `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\dcgpy `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\dcgpy `
    -DDCGP_BUILD_DCGP=no `
    -DDCGP_BUILD_DCGPY=yes ..

cmake --build . --target install --config Release
