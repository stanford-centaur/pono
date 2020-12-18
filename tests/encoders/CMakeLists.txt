include_directories("${PROJECT_SOURCE_DIR}/tests/encoders")

pono_add_test(test_btor2)
pono_add_test(test_coreir)
pono_add_test(test_smv)

# this is a long-running test
pono_add_test(test_btor2_ts_copy_equal)
