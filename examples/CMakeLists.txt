MACRO(ADD_EXAMPLE arg1)
	ADD_EXECUTABLE(${arg1} ${arg1}.cpp)
	TARGET_LINK_LIBRARIES(${arg1} dcgp Boost::program_options)
    target_compile_options(${arg1} PRIVATE "$<$<CONFIG:DEBUG>:${DCGP_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${DCGP_CXX_FLAGS_RELEASE}>")
    set_property(TARGET ${arg1} PROPERTY CXX_STANDARD 17)
    set_property(TARGET ${arg1} PROPERTY CXX_STANDARD_REQUIRED YES)
ENDMACRO(ADD_EXAMPLE)

# Symbolic regression examples (for the web docs)
ADD_EXAMPLE(symbolic_regression_1)
ADD_EXAMPLE(symbolic_regression_1mo)
ADD_EXAMPLE(symbolic_regression_2)
ADD_EXAMPLE(symbolic_regression_2mo)
ADD_EXAMPLE(symbolic_regression_3)
ADD_EXAMPLE(symbolic_regression_4)

# UDAs
ADD_EXAMPLE(es4cgp)
ADD_EXAMPLE(moes4cgp)
ADD_EXAMPLE(mes4cgp)
ADD_EXAMPLE(momes4cgp)


# Solving differential equations examples
ADD_EXAMPLE(tsoulos_ode1)
ADD_EXAMPLE(tsoulos_nlode3)

# Finding first integrals examples
ADD_EXAMPLE(hamiltonian_spring_mass)
ADD_EXAMPLE(hamiltonian_spring_mass_lipson)
ADD_EXAMPLE(hamiltonian_kepler)
