# Set Up Google Tests

# build testing infrastructure
# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

if (WITH_MSAT)
  add_definitions(-DWITH_MSAT)
endif()

if (WITH_CVC4)
  add_definitions(-DWITH_CVC4)
endif()

add_library(available-solvers SHARED
  "${PROJECT_SOURCE_DIR}/tests/available_solvers.cpp")
target_include_directories(available-solvers PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include")

add_executable(test_fts "${PROJECT_SOURCE_DIR}/tests/test_fts.cpp")
target_link_libraries(test_fts gtest_main)
target_link_libraries(test_fts cosa2)
target_link_libraries(test_fts available-solvers)
add_test(NAME test_fts COMMAND test_fts)
