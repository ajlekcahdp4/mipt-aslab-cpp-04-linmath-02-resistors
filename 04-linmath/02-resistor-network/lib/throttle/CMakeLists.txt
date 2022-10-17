set(LIB_SOURCES
  src/linear_solver.cc
)

add_library(throttle ${LIB_SOURCES})
target_include_directories(throttle PUBLIC include)
target_link_libraries(throttle PUBLIC range-v3::range-v3)

set(UNIT_TEST_SOURCES
  test/test_vector.cc
  test/test_contiguous_matrix.cc
  test/test_matrix.cc
  test/main.cc
)

if (ENABLE_GTEST)
  add_executable(unit_test ${UNIT_TEST_SOURCES})
  target_include_directories(unit_test PRIVATE src include)
  target_link_libraries(unit_test throttle ${GTEST_BOTH_LIBRARIES})
  gtest_discover_tests(unit_test)
endif()