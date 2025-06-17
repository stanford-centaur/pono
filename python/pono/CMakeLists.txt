configure_file(__init__.py.in __init__.py @ONLY)
configure_file(pono_imp.pxd pono_imp.pxd COPYONLY)

# Set Cython compilation flags.
set(CYTHON_FLAGS -3 --cplus --module-name pono._pono_impl)
if (WITH_COREIR OR WITH_COREIR_EXTERN)
  list(APPEND CYTHON_FLAGS --compile-time-env WITH_COREIR=ON)
else()
  list(APPEND CYTHON_FLAGS --compile-time-env WITH_COREIR=OFF)
endif()
if (WITH_MSAT_IC3IA)
  list(APPEND CYTHON_FLAGS --compile-time-env WITH_MSAT_IC3IA=ON)
else()
  list(APPEND CYTHON_FLAGS --compile-time-env WITH_MSAT_IC3IA=OFF)
endif()

add_custom_command(
  OUTPUT pono.cxx
  COMMAND
    "${Python_EXECUTABLE}" -m cython "${CMAKE_CURRENT_SOURCE_DIR}/pono.pyx" -o
    "${CMAKE_CURRENT_BINARY_DIR}/pono.cxx" ${CYTHON_FLAGS}
  DEPENDS pono.pyx
  COMMENT "Compiling Cython file ${CMAKE_CURRENT_SOURCE_DIR}/pono.pyx to ${CMAKE_CURRENT_BINARY_DIR}/pono.cxx"
)

Python_add_library(_pono_impl MODULE WITH_SOABI pono.cxx)
target_link_libraries(_pono_impl PUBLIC pono-lib)
# This has to be included because Cython uses the definitions from the smt-switch pxd files.
target_include_directories(_pono_impl PRIVATE "${SMT_SWITCH_DIR}/local/include/smt-switch")
