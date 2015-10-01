PROJECT(dCGP)

# We enable cmake test capability
ENABLE_TESTING()

# We set the minimum required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7)

# Allow CMAKE to find custom FindSOFTWARE.cmake fikes
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

# Set default build type to "Release", change it in the GUI if you need to build with debug.
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
		"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#Some screen output
MESSAGE(STATUS "OS detected: ${CMAKE_SYSTEM_NAME}")
MESSAGE(STATUS "CXX Compiler detected: ${CMAKE_CXX_COMPILER_ID}")

# Provides build options to CMake
# Build Option: when active the file main.cpp is built and linked to the PaGMO static library
OPTION(BUILD_MAIN "Build 'main.cpp'." ON)

# Build option: enable test set.
OPTION(BUILD_TESTS "Build test set." OFF)

# Build option: build examples
OPTION(BUILD_EXAMPLES "Build examples." OFF)

# Setting OS specific cmake variables
# 1 - OSX
#IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#	SET(CMAKE_MACOSX_RPATH OFF)
#	SET(CMAKE_OSX_ROOT "")
#	SET(CMAKE_OSX_DEPLOYMENT_TARGET "")
#	# Homebrew directory added to path
#	SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Initialise (empty) list of libraries to link against.
SET(MANDATORY_LIBRARIES "")



# Some screen output
MESSAGE(STATUS "CMake additional search path for libraries: ${CMAKE_LIBRARY_PATH}")

# We set compiler specific flags
INCLUDE(cmake_modules/set_compiler_flags.cmake)

# -----------------------------------------------------------------------------------
# We start finding third party dependencies

# 0 - Initial thread setup.
FIND_PACKAGE(Threads REQUIRED)
MESSAGE(STATUS "Thread library: ${CMAKE_THREAD_LIBS_INIT}")
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# 1 - Boost Setup.
IF(BUILD_TESTS)
	MESSAGE(STATUS "Linking unit tests to Boost.Test.")
	SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} unit_test_framework)
	IF(CMAKE_BUILD_TYPE STREQUAL "Release")
		MESSAGE(STATUS "Linking performance tests to Boost Timer and Boost.Chrono.")
		SET(REQUIRED_BOOST_LIBS ${REQUIRED_BOOST_LIBS} timer chrono system)
	ENDIF()
ENDIF()

MESSAGE(STATUS "Required Boost libraries: ${REQUIRED_BOOST_LIBS}")
# Timers are available since 1.48.0.
FIND_PACKAGE(Boost 1.48.0 REQUIRED COMPONENTS "${REQUIRED_BOOST_LIBS}")
MESSAGE(STATUS "Detected Boost version: ${Boost_VERSION}")
# Include system Boost headers.
MESSAGE(STATUS "Boost include dirs: ${Boost_INCLUDE_DIRS}")
MESSAGE(STATUS "Boost libraries: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# 2 - Piranha (sparse polynomial multiplication) Header only
FIND_PACKAGE(Piranha REQUIRED)
MESSAGE(STATUS "Piranha library found.")
MESSAGE(STATUS "Piranha include dir is: ${Piranha_INCLUDE_DIR}")
INCLUDE_DIRECTORIES(${Piranha_INCLUDE_DIR})

# 3 - Audi (automated differentiation) Header only
FIND_PACKAGE(Audi REQUIRED)
MESSAGE(STATUS "Audi library found.")
MESSAGE(STATUS "Audi include dir is: ${Piranha_INCLUDE_DIR}")
INCLUDE_DIRECTORIES(${Audi_INCLUDE_DIR})

# 2 - GMP setup (used by piranha)
FIND_PACKAGE(GMP REQUIRED)
MESSAGE(STATUS "GMP library found.")
MESSAGE(STATUS "GMP include dir is: ${GMP_INCLUDE_DIR}")
MESSAGE(STATUS "GMP library is: ${GMP_LIBRARIES}")
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})

# 3 - MPFR setup (used by piranha)
FIND_PACKAGE(MPFR REQUIRED)
MESSAGE(STATUS "MPFR library found.")
MESSAGE(STATUS "MPFR include dir is: ${MPFR_INCLUDE_DIR}")
MESSAGE(STATUS "MPFR library is: ${MPFR_LIBRARIES}")
MESSAGE(STATUS "Checking MPFR version.")
TRY_COMPILE(MPFR_VERSION_CHECK ${CMAKE_BINARY_DIR} "${CMAKE_SOURCE_DIR}/cmake_modules/mpfr_check_version.cpp"
	CMAKE_FLAGS "-DINCLUDE_DIRECTORIES:STRING=${MPFR_INCLUDE_DIR};${GMP_INCLUDE_DIR}")
IF(NOT MPFR_VERSION_CHECK)
	MESSAGE(FATAL_ERROR "Unsupported MPFR version, please upgrade.")
ENDIF(NOT MPFR_VERSION_CHECK)
MESSAGE(STATUS "MPFR version is ok.")
INCLUDE_DIRECTORIES(${MPFR_INCLUDE_DIR})

# NOTE: MPFR should be linked in before GMP, that's why we link GMP here.
SET(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

# Here we define the file list and make the static and dynamic library
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

# Build main and link it to static library.
IF(BUILD_MAIN)
	ADD_EXECUTABLE(main main.cpp)
	TARGET_LINK_LIBRARIES(main ${MANDATORY_LIBRARIES})
ENDIF(BUILD_MAIN)

# Build perf test in release and tests in debug
IF(BUILD_TESTS)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/tests")
ENDIF(BUILD_TESTS)

# Build pexamples
IF(BUILD_EXAMPLES)
	ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/examples")
ENDIF(BUILD_EXAMPLES)



