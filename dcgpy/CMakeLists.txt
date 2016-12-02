YACMA_PYTHON_MODULE(_core SHARED
	core.cpp docstrings.cpp
	# Python files.
	__init__.py
    test.py
    _graphviz_visualize.py
    _sympy_simplify.py
)

set_property(TARGET _core PROPERTY CXX_STANDARD 14)
set_property(TARGET _core PROPERTY CXX_STANDARD_REQUIRED YES)

# Setup the installation path.
set(DCGPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/dcgpy")

# Module needs the Python libs, plus the mandatory libraries.
if(WIN32)
    message(STATUS "Windows detected, linking also to the python library")
    set(MANDATORY_LIBRARIES ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

# Module needs the Python libs, plus the mandatory libraries.
TARGET_LINK_LIBRARIES(_core ${MANDATORY_LIBRARIES} ${Boost_SERIALIZATION_LIBRARY} ${Boost_PYTHON_LIBRARY})

INSTALL(TARGETS _core
	RUNTIME DESTINATION ${DCGPY_INSTALL_PATH}
	LIBRARY DESTINATION ${DCGPY_INSTALL_PATH}
)

# Add the global init file and the additional submodules.
set(DCGPY_PYTHON_MODULE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
    ${CMAKE_CURRENT_SOURCE_DIR}/test.py
    ${CMAKE_CURRENT_SOURCE_DIR}/_graphviz_visualize.py
    ${CMAKE_CURRENT_SOURCE_DIR}/_sympy_simplify.py
)

INSTALL(FILES ${DCGPY_PYTHON_MODULE_FILES} DESTINATION ${DCGPY_INSTALL_PATH})
