# Powershell script

# Install conda environment
conda config --set always_yes yes
conda create --name dcgpy cmake boost-cpp eigen pagmo-devel tbb audi symengine obake-devel cxx-compiler pybind11 python=3.8 sympy matplotlib
conda activate dcgpy

#build and run the dcgp tests
mkdir build
cd build
cmake `
    -G "Visual Studio 16 2019" `
    -A x64 `
    -DCMAKE_PREFIX_PATH=C:\Miniconda\envs\dcgpy `
    -DBoost_NO_BOOST_CMAKE=ON `
    -DCMAKE_BUILD_TYPE=Release `
    -DCMAKE_INSTALL_PREFIX=C:\Miniconda\envs\dcgpy `
    -DDCGP_BUILD_TESTS=no `
    -DDCGP_BUILD_EXAMPLES=no ..

cmake --build . --target install --config Release
