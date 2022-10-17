set(LIB_SOURCES
  src/resistor_network.cc
)

add_library(circuits ${LIB_SOURCES})
target_include_directories(circuits PUBLIC include)

set(UNIT_TEST_SOURCES
  test/main.cc
)

if (ENABLE_GTEST)
  add_executable(unit_test_circuits ${UNIT_TEST_SOURCES})
  target_include_directories(unit_test_circuits PRIVATE src include)
  target_link_libraries(unit_test_circuits throttle circuits ${GTEST_BOTH_LIBRARIES})
  gtest_discover_tests(unit_test_circuits)
endif()

target_link_libraries(circuits PUBLIC range-v3::range-v3 throttle)