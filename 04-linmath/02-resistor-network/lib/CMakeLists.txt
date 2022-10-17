add_subdirectory(throttle)
add_subdirectory(circuits)

include(FetchContent)
find_package(range-v3)

if (NOT range-v3_FOUND)
  FetchContent_Declare(
    range-v3
    GIT_REPOSITORY https://github.com/ericniebler/range-v3
    GIT_TAG master
  )

  if (NOT range-v3_POPULATED)
    FetchContent_Populate(range-v3)
    add_subdirectory(${range-v3_SOURCE_DIR} ${range-v3_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()

  FetchContent_MakeAvailable(range-v3)
endif()
