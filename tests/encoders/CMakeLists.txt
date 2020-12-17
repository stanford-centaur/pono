macro(pono_add_encoder_test name)
  add_executable(${name} "${PROJECT_SOURCE_DIR}/tests/encoders/${name}.cpp")
  target_link_libraries(${name} gtest_main)
  target_link_libraries(${name} pono-lib)
  add_test(NAME ${name} COMMAND ${name})
endmacro()

include_directories("${PROJECT_SOURCE_DIR}/tests/encoders")

pono_add_encoder_test(test_btor2)
pono_add_encoder_test(test_coreir)
pono_add_encoder_test(test_smv)

# this is a long-running test
pono_add_encoder_test(test_btor2_ts_copy_equal)
