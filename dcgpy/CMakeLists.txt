YACMA_PYTHON_MODULE(_core SHARED
	core.cpp docstrings.cpp
	# Python files.
	__init__.py
)

# Setup the installation path.
set(DCGPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/dcgpy")

# Module needs the Python libs, plus the mandatory libraries.
TARGET_LINK_LIBRARIES(_core ${MANDATORY_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS _core
	RUNTIME DESTINATION ${DCGPY_INSTALL_PATH}
	LIBRARY DESTINATION ${DCGPY_INSTALL_PATH}
)

# Add the global init file and the additional submodules.
INSTALL(FILES __init__.py DESTINATION ${DCGPY_INSTALL_PATH})
