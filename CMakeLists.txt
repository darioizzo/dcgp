# Main dcgp/dcgpy project version.
set(DCGP_PROJECT_VERSION 1.2.1)

# CMake version check.
cmake_minimum_required(VERSION 3.2)

# Module path setup.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules" "${CMAKE_SOURCE_DIR}/cmake_modules/yacma")

message(STATUS "System name: ${CMAKE_SYSTEM_NAME}")

# Set default build type to "Release".
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()

# Main build options: build dcgp or dcgpy. They cannot be on at the same time,
# and only one must be chosen.
option(DCGP_BUILD_DCGP "Build dcgp." ON)
option(DCGP_BUILD_DCGPY "Build dcgpy." OFF)

# Check consistency.
if(DCGP_BUILD_DCGP AND DCGP_BUILD_DCGPY)
    message(FATAL_ERROR "Please select whether to build dcgp or dcgpy: you cannot build them both at the same time.")
endif()

if((NOT DCGP_BUILD_DCGP) AND (NOT DCGP_BUILD_DCGPY))
    message(FATAL_ERROR "Please select if you want to build dcgp or dcgpy.")
endif()

if(DCGP_BUILD_DCGP)
    # Initial setup of a dcgp build.
    project(dcgp VERSION ${DCGP_PROJECT_VERSION})
    enable_testing()

    # Build option: enable test set.
    option(DCGP_BUILD_TESTS "Build test set." ON)

    # Build option: enable examples
    option(DCGP_BUILD_EXAMPLES "Build test set." OFF)

    # Build Option: when active the file main.cpp is built.
    option(DCGP_BUILD_MAIN "Build 'main.cpp'." OFF)
else()
    # Initial setup of a dcgpy build.
    project(dcgpy VERSION ${DCGP_PROJECT_VERSION})
endif()

# Common general bits.

# Initial setup of compiler flags.
include(YACMACompilerLinkerSettings)

# Threading setup.
include(YACMAThreadingSetup)

# Assemble the flags.
set(DCGP_CXX_FLAGS_DEBUG ${YACMA_CXX_FLAGS} ${YACMA_CXX_FLAGS_DEBUG} ${YACMA_THREADING_CXX_FLAGS})
set(DCGP_CXX_FLAGS_RELEASE ${YACMA_CXX_FLAGS} ${YACMA_THREADING_CXX_FLAGS})
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND YACMA_COMPILER_IS_CLANGXX)
  message(STATUS "Clang compiler on OSX detected, setting the standard library to 'libc++'.")
  list(APPEND DCGP_CXX_FLAGS_DEBUG "-stdlib=libc++")
  list(APPEND DCGP_CXX_FLAGS_RELEASE "-stdlib=libc++")
endif()
if(YACMA_COMPILER_IS_MSVC)
  # Disable the idiotic minmax macros on MSVC, some annoying warnings,
  # and enable the bigobj option.
  list(APPEND DCGP_CXX_FLAGS_DEBUG "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj")
  list(APPEND DCGP_CXX_FLAGS_RELEASE "-DNOMINMAX" "/wd4459" "/wd4127" "/wd4702" "/bigobj")
endif()
if(YACMA_COMPILER_IS_INTELXX)
  # NOTE: on MSVC we use the push/pop pragmas, but they do not seem to work on Intel (the pragmas
  # in icc influence the behaviour at instantiation point, not at definition point).
  # These warnings are useful in principle, but they are generated a lot from cereal and we have no
  # way of disabling them selectively. Just rely on the other compilers to provde good diagnostic.
  list(APPEND DCGP_CXX_FLAGS_DEBUG "-diag-disable" "2259,1682,68")
  list(APPEND DCGP_CXX_FLAGS_RELEASE "-diag-disable" "2259,1682,68")
endif()
if(MINGW)
	# Flag needed to deal with big binaries in MinGW.
	message(STATUS "Enabling the '-Wa,-mbig-obj' flag in MinGW builds.")
	list(APPEND DCGP_CXX_FLAGS_DEBUG "-Wa,-mbig-obj")
	list(APPEND DCGP_CXX_FLAGS_RELEASE "-Wa,-mbig-obj")
endif()

# Some flags that generate warnings due to Eigen or Piranha obsolete versions
list(REMOVE_ITEM DCGP_CXX_FLAGS_DEBUG "-Wduplicated-branches")
list(REMOVE_ITEM DCGP_CXX_FLAGS_DEBUG "-Wold-style-cast")

# Creates the file config.hpp
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/include/dcgp/config.hpp" @ONLY)

if(DCGP_BUILD_DCGPY)
    include(YACMAPythonSetup)
    # Python version check.
    if(${PYTHON_VERSION_MAJOR} LESS 2 OR (${PYTHON_VERSION_MAJOR} EQUAL 2 AND ${PYTHON_VERSION_MINOR} LESS 7))
        message(FATAL_ERROR "Minimum supported Python version is 2.7.")
    endif()

    # NOTE: for the time being, require that dcgp/dcgpy versions are matching exactly.
    find_package(dcgp ${DCGP_PROJECT_VERSION} EXACT REQUIRED)

    # Numpy headers are needed to convert to/from std container
    find_package(NumPy REQUIRED)

    # Internal variable that will be used to tell DCGPFindBoost to locate Boost.Python.
    set(_DCGP_FIND_BOOST_PYTHON TRUE)
endif()

# Boost setup (common to dcgp/dcgpy).
include(DCGPFindBoost)

if(DCGP_BUILD_DCGP)
    # Eigen setup
    find_package(Eigen3 REQUIRED)
    message(STATUS "Eigen include directory: ${EIGEN3_INCLUDE_DIR}")
    message(STATUS "Eigen version detected: ${EIGEN3_VERSION}")
    include_directories(${EIGEN3_INCLUDE_DIR})

    # Piranha setup.
    find_package(Piranha REQUIRED)
    message(STATUS "Piranha library found.")
    message(STATUS "Piranha include dir is: ${Piranha_INCLUDE_DIR}")
    include_directories(${Piranha_INCLUDE_DIR})

    # GMP setup.
    find_package(GMP REQUIRED)
    message(STATUS "GMP library found.")
    message(STATUS "GMP include dir is: ${GMP_INCLUDE_DIR}")
    message(STATUS "GMP library is: ${GMP_LIBRARIES}")

    # MPFR setup.
    find_package(MPFR REQUIRED)
    message(STATUS "MPFR library found.")
    message(STATUS "MPFR include dir is: ${MPFR_INCLUDE_DIR}")
    message(STATUS "MPFR library is: ${MPFR_LIBRARIES}")
    message(STATUS "Checking MPFR version.")
    try_compile(MPFR_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/mpfr_check_version.cpp"
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MPFR_INCLUDE_DIR};${GMP_INCLUDE_DIR}")
    if(NOT MPFR_VERSION_CHECK)
        message(FATAL_ERROR "Unsupported MPFR version, please upgrade.")
    endif()
    message(STATUS "MPFR version is ok.")

    # Audi setup
    find_package(Audi REQUIRED)
    message(STATUS "AUDI include dir is: ${AUDI_INCLUDE_DIRS}")

    # Pagmo setup
    find_package(Pagmo REQUIRED)
    message(STATUS "Pagmo dir is: ${Pagmo_DIR}")

  
    # TBB setup
    find_package(TBB REQUIRED)
    message(STATUS "TBB include dir is: ${TBB_INCLUDE_DIRS}")
    message(STATUS "TBB library is: ${TBB_LIBRARIES}")


    # Setup of the header-only dcgp library.
    add_library(dcgp INTERFACE)
    target_link_libraries(dcgp INTERFACE Threads::Threads Boost::boost Boost::serialization Eigen3::eigen3 MPFR::MPFR GMP::GMP TBB::tbb)
    target_link_libraries(dcgp INTERFACE Audi::audi Pagmo::pagmo)


    # This sets up the include directory to be different if we build
    target_include_directories(dcgp INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)

    # These lines are temporary with piranha v0.10, may be simplified when this is updated
    target_include_directories(dcgp INTERFACE ${Piranha_INCLUDE_DIR})

    # Build main
    if(DCGP_BUILD_MAIN)
        add_executable(main main.cpp)
        target_link_libraries(main dcgp)
        set_property(TARGET main PROPERTY CXX_STANDARD 14)
        set_property(TARGET main PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET main PROPERTY CXX_EXTENSIONS NO)
    endif()

    # Builds the tests (performance in Release, the others in Debug)
    if(DCGP_BUILD_TESTS)
        add_subdirectory("${CMAKE_SOURCE_DIR}/tests")
    endif()

    # Builds the examples
    if(DCGP_BUILD_EXAMPLES)
        add_subdirectory("${CMAKE_SOURCE_DIR}/examples")
    endif()

    # Configure the doc files.
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/doc/sphinx/conf.py" @ONLY)

    # Setup of the export.
    install(TARGETS dcgp EXPORT dcgp_export)

    # Setup cmake config files
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dcgp-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/dcgp-config.cmake" @ONLY)

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dcgp-config.cmake" DESTINATION "lib/cmake/dcgp")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindGMP.cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindMPFR.cmake"
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindPiranha.cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/DCGPFindBoost.cmake"
            "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindAudi.cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindEigen3.cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/FindTBB.cmake" DESTINATION "lib/cmake/dcgp")
    install(EXPORT dcgp_export NAMESPACE Dcgp:: DESTINATION lib/cmake/dcgp)
    # Take care of versioning.
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/dcgp-config-version.cmake" VERSION ${dcgp_VERSION}
            COMPATIBILITY ExactVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/dcgp-config-version.cmake" DESTINATION "lib/cmake/dcgp")

    # Do the actual library installation.
    install(DIRECTORY include/ DESTINATION include)

endif()

if(DCGP_BUILD_DCGPY)
    # Configure version
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/dcgpy/_version.py.in" "${CMAKE_CURRENT_SOURCE_DIR}/dcgpy/_version.py" @ONLY)
    # Add the dcgpy subdirectory.
    add_subdirectory("${CMAKE_SOURCE_DIR}/dcgpy")
    if(MINGW)
        message(STATUS "Creating the files for the generation of a binary wheel.")
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py" @ONLY)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
    endif()
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # NOTE: this is necessary on linux but harmful on mingw.
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.cfg" @ONLY)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/wheel/setup.py" @ONLY)
    endif()
endif()

