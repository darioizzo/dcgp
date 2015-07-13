ADD_EXECUTABLE(test_compute_f test_compute_f.cpp)
TARGET_LINK_LIBRARIES(test_compute_f ${MANDATORY_LIBRARIES} dcgp_s)
ADD_TEST(test_compute_f test_compute_f)