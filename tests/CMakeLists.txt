# Add the unit test library.
# NOTE: probably we should have a check in the build system to determine
# whether Boost libraries are static or dynamic.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)

MACRO(ADD_DCGP_TESTCASE arg1)
    IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
        ADD_EXECUTABLE(${arg1} ${arg1}.cpp)
        TARGET_LINK_LIBRARIES(${arg1} dcgp Boost::unit_test_framework)
        target_compile_options(${arg1} PRIVATE "$<$<CONFIG:DEBUG>:${DCGP_CXX_FLAGS_DEBUG}>")		
        set_property(TARGET ${arg1} PROPERTY CXX_STANDARD 14)
        set_property(TARGET ${arg1} PROPERTY CXX_STANDARD_REQUIRED YES)
        ADD_TEST(${arg1} ${arg1})
    ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
ENDMACRO(ADD_DCGP_TESTCASE)

MACRO(ADD_DCGP_PERFORMANCE_TESTCASE arg1)
    IF(CMAKE_BUILD_TYPE STREQUAL "Release")
        ADD_EXECUTABLE("${arg1}_perf" "${arg1}_perf.cpp")
        TARGET_LINK_LIBRARIES("${arg1}_perf" dcgp Boost::timer Boost::chrono Boost::system Boost::unit_test_framework)
        target_compile_options("${arg1}_perf" PRIVATE "$<$<CONFIG:RELEASE>:${DCGP_CXX_FLAGS_RELEASE}>")		
        set_property(TARGET "${arg1}_perf" PROPERTY CXX_STANDARD 14)
        set_property(TARGET "${arg1}_perf" PROPERTY CXX_STANDARD_REQUIRED YES)
        ADD_TEST("${arg1}_perf" "${arg1}_perf")
    ENDIF(CMAKE_BUILD_TYPE STREQUAL "Release")
ENDMACRO(ADD_DCGP_PERFORMANCE_TESTCASE)

ADD_DCGP_TESTCASE(expression)
ADD_DCGP_TESTCASE(differentiate)
ADD_DCGP_TESTCASE(expression_ann)
ADD_DCGP_TESTCASE(wrapped_functions)
ADD_DCGP_TESTCASE(rng)


ADD_DCGP_PERFORMANCE_TESTCASE(function_calls)
ADD_DCGP_PERFORMANCE_TESTCASE(compute)
ADD_DCGP_PERFORMANCE_TESTCASE(loss)
ADD_DCGP_PERFORMANCE_TESTCASE(mutate)
ADD_DCGP_PERFORMANCE_TESTCASE(differentiate)
ADD_DCGP_PERFORMANCE_TESTCASE(expression_ann)


