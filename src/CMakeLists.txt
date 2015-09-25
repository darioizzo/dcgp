# Keplerian_toolbox lib source files.
SET(dCGP_LIB_SRC_LIST
	${CMAKE_CURRENT_SOURCE_DIR}/expression.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/function_set.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/rng.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/wrapped_functions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fitness_functions.cpp
)

#Build Static Library
ADD_LIBRARY(dcgp_s STATIC ${dCGP_LIB_SRC_LIST})

