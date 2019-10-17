MACRO(ADD_EXAMPLE arg1)
	ADD_EXECUTABLE(${arg1} ${arg1}.cpp)
	TARGET_LINK_LIBRARIES(${arg1} dcgp)
    target_compile_options(${arg1} PUBLIC "$<$<CONFIG:DEBUG>:${YACMA_CXX_FLAGS_DEBUG}>")
    target_compile_options(${arg1} PUBLIC "$<$<CONFIG:RELEASE>:${YACMA_CXX_FLAGS}>")
    set_property(TARGET ${arg1} PROPERTY CXX_STANDARD 14)
    set_property(TARGET ${arg1} PROPERTY CXX_STANDARD_REQUIRED YES)
ENDMACRO(ADD_EXAMPLE)

# Symbolic regression examples
ADD_EXAMPLE(symbolic_regression_1)

# Solving differential equations examples
ADD_EXAMPLE(tsoulos_ode1)
ADD_EXAMPLE(tsoulos_nlode3)

# Finding first integrals examples
ADD_EXAMPLE(hamiltonian_spring_mass)
ADD_EXAMPLE(hamiltonian_spring_mass_lipson)
ADD_EXAMPLE(hamiltonian_kepler)
