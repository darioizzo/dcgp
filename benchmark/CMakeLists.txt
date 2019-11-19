function(ADD_DCGP_BENCHMARK arg1)
  add_executable(${arg1} ${arg1}.cpp)
  target_link_libraries(${arg1} dcgp Boost::timer Boost::chrono Boost::system)
  target_compile_options(${arg1} PRIVATE
    "$<$<CONFIG:Debug>:${DCGP_CXX_FLAGS_DEBUG}>"
    "$<$<CONFIG:Release>:${DCGP_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:RelWithDebInfo>:${DCGP_CXX_FLAGS_RELEASE}>"
    "$<$<CONFIG:MinSizeRel>:${DCGP_CXX_FLAGS_RELEASE}>"
  )
  set_property(TARGET ${arg1} PROPERTY CXX_STANDARD 17)
  set_property(TARGET ${arg1} PROPERTY CXX_STANDARD_REQUIRED YES)
  set_property(TARGET ${arg1} PROPERTY CXX_EXTENSIONS NO)
endfunction()

ADD_DCGP_BENCHMARK(function_calls)
ADD_DCGP_BENCHMARK(compute)
ADD_DCGP_BENCHMARK(loss)
ADD_DCGP_BENCHMARK(mutate)
ADD_DCGP_BENCHMARK(differentiate)
ADD_DCGP_BENCHMARK(expression_ann)
ADD_DCGP_BENCHMARK(symbolic_regression)


