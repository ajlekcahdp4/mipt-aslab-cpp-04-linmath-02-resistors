set(NET_SOURCES
  src/network.cc
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
target_link_libraries(network throttle circuits Boost::program_options scn)

install(TARGETS network DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

set(COMP_SOURCES
  src/roughly_compare.cc
)


add_executable(comp ${COMP_SOURCES})
target_link_libraries(comp PRIVATE throttle Boost::program_options)
target_link_libraries(comp PRIVATE scn)
install(TARGETS comp DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

if(BASH_PROGRAM)
    add_test(NAME test.network COMMAND ${BASH_PROGRAM} ${CMAKE_CURRENT_SOURCE_DIR}/test.sh "$<TARGET_FILE:network>" ${CMAKE_CURRENT_SOURCE_DIR} "$<TARGET_FILE:comp>")
endif()

