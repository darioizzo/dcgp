IF(AUDI_INCLUDE_DIR)
	# Already in cache, be silent
	SET(AUDI_FIND_QUIETLY TRUE)
ENDIF(AUDI_INCLUDE_DIR)

FIND_PATH(AUDI_INCLUDE_DIR NAMES audi/audi.hpp)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Audi DEFAULT_MSG AUDI_INCLUDE_DIR)

MARK_AS_ADVANCED(AUDI_INCLUDE_DIR)

# Added July 2018. Without appveyor (only) would complain for a missing target Audi::audi. Not clear why Piranha does not have it ... maybe it should
if(AUDI_INCLUDE_DIR AND NOT TARGET Audi::audi)
    add_library(Audi::audi INTERFACE IMPORTED)
    set_target_properties(Audi::audi PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${AUDI_INCLUDE_DIR}")
endif()
