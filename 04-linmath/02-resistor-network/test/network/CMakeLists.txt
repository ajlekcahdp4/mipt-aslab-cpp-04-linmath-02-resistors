set(NET_SOURCES
  src/network.cc
  src/resistor_network.cc
)

FetchContent_Declare(
  scn
  GIT_REPOSITORY https://github.com/eliaskosunen/scnlib
  GIT_TAG master
)

if (NOT scn_POPULATED)
  FetchContent_Populate(scn)
  add_subdirectory(${scn_SOURCE_DIR} ${scn_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

add_executable(network ${NET_SOURCES})
target_link_libraries(network throttle Boost::program_options scn)

install(TARGETS network DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# if(BASH_PROGRAM)
#   add_test(NAME test.determinant COMMAND ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/test.sh "$<TARGET_FILE:determinant>" ${CMAKE_CURRENT_SOURCE_DIR} "$<TARGET_FILE:comp>")
# endif()

