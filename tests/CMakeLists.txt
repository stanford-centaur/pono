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

add_library(pono-test-lib "${PONO_LIB_TYPE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/common_ts.cpp"
)

# INCLUDE_DIRS set in top-level CMakeLists.txt
target_include_directories(pono-test-lib PUBLIC "${INCLUDE_DIRS}")
target_link_libraries(pono-test-lib pono-lib)

macro(pono_add_test name)
  add_executable(${name} "${CMAKE_CURRENT_SOURCE_DIR}/${name}.cpp")
  target_link_libraries(${name} gtest_main)
  target_link_libraries(${name} pono-test-lib) # also includes pono-lib because it's linked
  add_test(NAME ${name} COMMAND ${name})
endmacro()

include_directories("${PROJECT_SOURCE_DIR}/tests/")

pono_add_test(test_ts)
pono_add_test(test_unroller)
pono_add_test(test_modifiers)
pono_add_test(test_engines)
pono_add_test(test_utils)
pono_add_test(test_uf)
pono_add_test(test_witness)
pono_add_test(test_control_signal_modifiers)
pono_add_test(test_coi)
pono_add_test(test_ts_replace_terms)
pono_add_test(test_ic3)
pono_add_test(test_ic3bits)
pono_add_test(test_ic3ia)
pono_add_test(test_ic3sa)
pono_add_test(test_msat_ic3ia)
pono_add_test(test_ceg_prophecy_arrays)
pono_add_test(test_cegar_ops_uf)
pono_add_test(test_cegar_values)
pono_add_test(test_term_analysis)
pono_add_test(test_walkers)
pono_add_test(test_pseudo_init_and_prop)
pono_add_test(test_promote_inputvars)
pono_add_test(test_partial_model)

add_subdirectory(encoders)
