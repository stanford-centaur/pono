if(POLICY CMP0057)
  # For cmake >= 3.3 this policy changed the behavior of IN_LIST
  # if the policy exists, we use the NEW behavior
  cmake_policy(SET CMP0057 NEW)
endif()

find_package(Cython 0.29 REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/deps/smt-switch/python")

add_cython_target(cosa2 CXX)

add_library(cosa2 MODULE ${cosa2})

target_link_libraries(cosa2 cosa2-lib)
target_include_directories(cosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(cosa2)

configure_file(setup.py.in setup.py)

install(TARGETS cosa2 DESTINATION lib)
