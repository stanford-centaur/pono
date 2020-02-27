if(POLICY CMP0057)
  # For cmake >= 3.3 this policy changed the behavior of IN_LIST
  # if the policy exists, we use the NEW behavior
  cmake_policy(SET CMP0057 NEW)
endif()

find_package(Cython 0.29 REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/deps/smt-switch/python")

add_cython_target(pycosa2 CXX)

add_library(pycosa2 MODULE ${pycosa2})

target_link_libraries(pycosa2 cosa2)
target_include_directories(pycosa2 PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(pycosa2)

install(TARGETS pycosa2 DESTINATION lib)
