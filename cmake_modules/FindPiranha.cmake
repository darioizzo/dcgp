IF(Piranha_INCLUDE_DIR)
	# Already in cache, be silent
	SET(Piranha_FIND_QUIETLY TRUE)
ENDIF(Piranha_INCLUDE_DIR)

FIND_PATH(Piranha_INCLUDE_DIR NAMES piranha/piranha.hpp)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Piranha DEFAULT_MSG Piranha_INCLUDE_DIR)

MARK_AS_ADVANCED(Piranha_INCLUDE_DIR)
