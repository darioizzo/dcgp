PROJECT(dCGP)

# We enable cmake test capability
ENABLE_TESTING()

# We set the minimum required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#Some screen output
MESSAGE(STATUS "OS detected: ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "CXX Compiler detected: ${CMAKE_CXX_COMPILER_ID}")

# Setting OS specific cmake variables
# 1 - OSX
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	SET(CMAKE_MACOSX_RPATH OFF)
	SET(CMAKE_OSX_ROOT "")
	SET(CMAKE_OSX_DEPLOYMENT_TARGET "")
	# Homebrew directory added to path
	SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Some screen output
MESSAGE(STATUS "CMake additional search path for libraries: ${CMAKE_LIBRARY_PATH}")

# We set compiler specific flags
INCLUDE(cmake_scripts/set_compiler_flags.cmake)

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build option: enable test set.
OPTION(ENABLE_TESTS "Build test set." OFF)

# Finding the boost libraries needed for the keplerian_toolbox
# SET(REQUIRED_BOOST_LIBS serialization)
# MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
# FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
# MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
# MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
# MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
# INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# Initialise (empty) list of libraries to link against.
SET(LIBRARIES_4_STATIC "")
SET(LIBRARIES_4_DYNAMIC "")

# Define the libraries to link against.
SET(LIBRARIES_4_STATIC ${LIBRARIES_4_STATIC})
SET(LIBRARIES_4_DYNAMIC ${LIBRARIES_4_DYNAMIC} ${LIBRARIES_4_STATIC})

# Here we define the file list and make the static and dynamic library
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Build main and link it to static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main dcgp_s ${LIBRARIES_4_STATIC})
ENDIF(BUILD_MAIN)

IF(ENABLE_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(ENABLE_TESTS)


