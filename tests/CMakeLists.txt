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

include_directories("${PROJECT_SOURCE_DIR}/tests/encoders")
include_directories("${PROJECT_SOURCE_DIR}/tests/")

add_executable(test_ts "${PROJECT_SOURCE_DIR}/tests/test_ts.cpp")
target_link_libraries(test_ts gtest_main)
target_link_libraries(test_ts pono-lib)
add_test(NAME test_ts COMMAND test_ts)

add_executable(test_unroller "${PROJECT_SOURCE_DIR}/tests/test_unroller.cpp")
target_link_libraries(test_unroller gtest_main)
target_link_libraries(test_unroller pono-lib)
add_test(NAME test_unroller COMMAND test_unroller)

add_executable(test_modifiers "${PROJECT_SOURCE_DIR}/tests/test_modifiers.cpp")
target_link_libraries(test_modifiers gtest_main)
target_link_libraries(test_modifiers pono-lib)
add_test(NAME test_modifiers COMMAND test_modifiers)

add_executable(test_engines "${PROJECT_SOURCE_DIR}/tests/test_engines.cpp")
target_link_libraries(test_engines gtest_main)
target_link_libraries(test_engines pono-lib)
add_test(NAME test_engines COMMAND test_engines)

add_executable(test_utils "${PROJECT_SOURCE_DIR}/tests/test_utils.cpp")
target_link_libraries(test_utils gtest_main)
target_link_libraries(test_utils pono-lib)
add_test(NAME test_utils COMMAND test_utils)

add_executable(test_uf "${PROJECT_SOURCE_DIR}/tests/test_uf.cpp")
target_link_libraries(test_uf gtest_main)
target_link_libraries(test_uf pono-lib)
add_test(NAME test_uf COMMAND test_uf)

add_executable(test_witness "${PROJECT_SOURCE_DIR}/tests/test_witness.cpp")
target_link_libraries(test_witness gtest_main)
target_link_libraries(test_witness pono-lib)
add_test(NAME test_witness COMMAND test_witness)

add_executable(test_control_signal_modifiers "${PROJECT_SOURCE_DIR}/tests/test_control_signal_modifiers.cpp")
target_link_libraries(test_control_signal_modifiers gtest_main)
target_link_libraries(test_control_signal_modifiers pono-lib)
add_test(NAME test_control_signal_modifiers COMMAND test_control_signal_modifiers)

add_executable(test_btor2_ts_copy_equal "${PROJECT_SOURCE_DIR}/tests/test_btor2_ts_copy_equal.cpp")
target_link_libraries(test_btor2_ts_copy_equal gtest_main)
target_link_libraries(test_btor2_ts_copy_equal pono-lib)
add_test(NAME test_btor2_ts_copy_equal COMMAND test_btor2_ts_copy_equal)

add_executable(test_coreir "${PROJECT_SOURCE_DIR}/tests/encoders/test_coreir.cpp")
target_link_libraries(test_coreir gtest_main)
target_link_libraries(test_coreir pono-lib)
add_test(NAME test_coreir COMMAND test_coreir)

add_executable(test_ic3 "${PROJECT_SOURCE_DIR}/tests/test_ic3.cpp")
target_link_libraries(test_ic3 gtest_main)
target_link_libraries(test_ic3 pono-lib)
add_test(NAME test_ic3 COMMAND test_ic3)

add_executable(test_btor2 "${PROJECT_SOURCE_DIR}/tests/encoders/test_btor2.cpp")
target_link_libraries(test_btor2 gtest_main)
target_link_libraries(test_btor2 pono-lib)
add_test(NAME test_btor2 COMMAND test_btor2)

add_executable(test_ceg_prophecy_arrays "${PROJECT_SOURCE_DIR}/tests/test_ceg_prophecy_arrays.cpp")
target_link_libraries(test_ceg_prophecy_arrays gtest_main)
target_link_libraries(test_ceg_prophecy_arrays pono-lib)
add_test(NAME test_ceg_prophecy_arrays COMMAND test_ceg_prophecy_arrays)
