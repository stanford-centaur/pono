if(POLICY CMP0057)
  # For cmake >= 3.3 this policy changed the behavior of IN_LIST
  # if the policy exists, we use the NEW behavior
  cmake_policy(SET CMP0057 NEW)
endif()

find_package(Cython 0.29 REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/deps/smt-switch/python")

if (WITH_COREIR)
  set(CYTHON_FLAGS "${CYTHON_FLAGS}" "--compile-time-env" "WITH_COREIR=ON")
else()
  set(CYTHON_FLAGS "${CYTHON_FLAGS}" "--compile-time-env" "WITH_COREIR=OFF")
endif()

add_cython_target(pono CXX)

add_library(pono MODULE ${pono})

target_link_libraries(pono pono-lib)
target_include_directories(pono PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(pono)

configure_file(setup.py.in setup.py)

install(TARGETS pono DESTINATION lib)
