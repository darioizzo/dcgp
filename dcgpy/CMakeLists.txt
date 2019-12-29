if(UNIX)
# On unix platforms, we set the rpath of the dcgpy libraries
# to the lib subdirectory of the install prefix. This allows
# to install dcgpy in nonstandard locations (e.g., ~/.local)
# without having to muck around with library search paths.
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

# Finding pygmo
find_package(Pygmo REQUIRED)

# The name of the BP target to link to depends on the Boost and Python versions.
# NOTE: since Boost 1.67, the naming of the Boost.Python library has changed to include the
# major and minor python version as a suffix. See the release notes:
# https://www.boost.org/users/history/version_1_67_0.html
if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 66))
    set(DCGPY_BP_TARGET "Boost::python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}")
else()
    if(${PYTHON_VERSION_MAJOR} EQUAL 2)
        set(DCGPY_BP_TARGET "Boost::python")
    else()
        set(DCGPY_BP_TARGET "Boost::python3")
    endif()
endif()

# Boost.Python version 1.64 (and possibly later?) has issues with default hidden visibility.
if(${Boost_MAJOR_VERSION} GREATER 1 OR (${Boost_MAJOR_VERSION} EQUAL 1 AND ${Boost_MINOR_VERSION} GREATER 63))
list(FIND DCGP_CXX_FLAGS_DEBUG "-fvisibility=hidden" _DCGPY_HAS_VHIDDEN)
if (NOT ${_DCGPY_HAS_VHIDDEN} EQUAL -1)
	message(STATUS "Removing the '-fvisibility=hidden' flag for Boost.Python > 1.63.")
	list(REMOVE_ITEM DCGP_CXX_FLAGS_DEBUG "-fvisibility=hidden")
	list(REMOVE_ITEM DCGP_CXX_FLAGS_RELEASE "-fvisibility=hidden")
endif()
endif()

# Setup of the header-only dcgpy library.
add_library(dcgpy INTERFACE)
target_link_libraries(dcgpy INTERFACE ${DCGPY_BP_TARGET} Boost::disable_autolinking YACMA::PythonIncludeDir NumPy::numpy Dcgp::dcgp Pygmo::pygmo)
target_include_directories(dcgpy INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../>
$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include/dcgpy> $<INSTALL_INTERFACE:include>)

# Extra flags to be appended in CI experiments
set(DCGP_CXX_FLAGS_EXTRA)
mark_as_advanced(DCGP_CXX_FLAGS_EXTRA)
list(APPEND DCGP_CXX_FLAGS_RELEASE ${DCGP_CXX_FLAGS_EXTRA})

# Setup of the dcgpy core module.
YACMA_PYTHON_MODULE(core 
    core.cpp 
    expose_expressions.cpp
    expose_kernels.cpp
    expose_symbolic_regression.cpp
    docstrings.cpp)
target_link_libraries(core PRIVATE ${DCGPY_BP_TARGET} dcgpy)
target_compile_options(core PRIVATE "$<$<CONFIG:DEBUG>:${DCGP_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${DCGP_CXX_FLAGS_RELEASE}>")
set_property(TARGET core PROPERTY CXX_STANDARD 17)
set_property(TARGET core PROPERTY CXX_STANDARD_REQUIRED YES)
set_property(TARGET core PROPERTY CXX_EXTENSIONS NO)

# Setup the installation path.
set(DCGPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/dcgpy")
install(TARGETS core
RUNTIME DESTINATION ${DCGPY_INSTALL_PATH}
LIBRARY DESTINATION ${DCGPY_INSTALL_PATH}
)

# Add the Python files.
install(FILES __init__.py _version.py test.py _visualize.py _simplify.py _encode.py DESTINATION ${DCGPY_INSTALL_PATH})
