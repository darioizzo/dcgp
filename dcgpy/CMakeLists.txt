# Python version check.
if(${PYTHON_VERSION_MAJOR} LESS 3)
    message(FATAL_ERROR "Minimum supported Python version is 3.0")
endif()

# Setting the rpath
if(UNIX)
    # On unix platforms, we set the rpath of the pygmo libraries
    # to the lib subdirectory of the install prefix. This allows
    # to install pygmo in nonstandard locations (e.g., ~/.local)
    # without having to muck around with library search paths.
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Setup of the dcgpy core module.
YACMA_PYTHON_MODULE(core 
    core.cpp 
    expose_expressions.cpp
    expose_kernels.cpp
    #expose_symbolic_regression.cpp
    docstrings.cpp
    )

# Finding dcgp 
find_package(dcgp ${DCGP_PROJECT_VERSION} EXACT REQUIRED)

# Finding pybind11
find_package(pybind11 REQUIRED)

# Setup of the header-only dcgpy library.
add_library(dcgpy INTERFACE)
target_link_libraries(dcgpy INTERFACE ${DCGPY_BP_TARGET} Boost::disable_autolinking Dcgp::dcgp)
target_include_directories(dcgpy INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/dcgpy> $<INSTALL_INTERFACE:include>)

# Extra flags to be appended in CI experiments
set(DCGP_CXX_FLAGS_EXTRA)
mark_as_advanced(DCGP_CXX_FLAGS_EXTRA)
list(APPEND DCGP_CXX_FLAGS_RELEASE ${DCGP_CXX_FLAGS_EXTRA})

# Setup of the dcgpy core module.
target_link_libraries(core PRIVATE ${DCGPY_BP_TARGET} dcgpy)
target_compile_options(core PRIVATE "$<$<CONFIG:DEBUG>:${DCGP_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${DCGP_CXX_FLAGS_RELEASE}>")
target_include_directories(core SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")

set_property(TARGET core PROPERTY CXX_STANDARD 17)
set_property(TARGET core PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(DCGPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/dcgpy")
install(TARGETS core
RUNTIME DESTINATION ${DCGPY_INSTALL_PATH}
LIBRARY DESTINATION ${DCGPY_INSTALL_PATH}
)


# Helper file with version number.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/_version.py" @ONLY)   

# Add the Python files.
install(FILES __init__.py _version.py test.py _visualize.py _simplify.py _encode.py DESTINATION ${DCGPY_INSTALL_PATH})
