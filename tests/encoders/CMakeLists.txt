include_directories("${PROJECT_SOURCE_DIR}/tests/encoders")

pono_add_test(test_btor2)
pono_add_test(test_coreir)
pono_add_test(test_smv)
pono_add_test(test_systemverilog_types)
pono_add_test(test_systemverilog_operators)
pono_add_test(test_systemverilog_statements)
pono_add_test(test_systemverilog_hierarchy)
pono_add_test(test_systemverilog_sva)
pono_add_test(test_systemverilog_unsupported)
pono_add_test(test_systemverilog_integration)

# this is a long-running test
pono_add_test(test_btor2_ts_copy_equal)
