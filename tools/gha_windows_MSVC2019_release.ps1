# Powershell script
# Install conda environment
conda config --set always_yes yes
conda create --name dcgpy cmake boost-cpp eigen pagmo-devel tbb audi symengine obake-devel cxx-compiler pybind11 clang ninja
conda activate dcgpy

# Define shell variables for clang
& "C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Auxiliary\Build\vcvars64.bat && set > %temp%\vcvars.txt"

Get-Content "$env:temp\vcvars.txt" | Foreach-Object {
  if ($_ -match "^(.*?)=(.*)$") {
    Set-Content "env:\$($matches[1])" $matches[2]
  }
}

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

cmake --build . --config Release
ctest -j4 -V -C Release