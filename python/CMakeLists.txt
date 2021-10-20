# relying on cmake files from scikit-build in smt-switch deps
if (NOT EXISTS "${SMT_SWITCH_DIR}/deps/scikit-build/skbuild/resources/cmake")
  message(FATAL_ERROR "Missing CMake files in smt-switch. Remember to set it up with --python")
else()
  list(APPEND CMAKE_MODULE_PATH ${SMT_SWITCH_DIR}/deps/scikit-build/skbuild/resources/cmake)
endif()

if(POLICY CMP0057)
  # For cmake >= 3.3 this policy changed the behavior of IN_LIST
  # if the policy exists, we use the NEW behavior
  cmake_policy(SET CMP0057 NEW)
endif()

# We currently use FindPythonInterp even though it is deprecated since 3.12
# This is because the scikit-build files still use this version and it will
# not interact well with the latest Python finding cmake modules
# https://cmake.org/cmake/help/v3.12/module/FindPython.html
# in the future, we can switch to FindPython3 once it has become more standard
# i.e. when the following issue is resolved:
# https://github.com/scikit-build/scikit-build/issues/506
find_package(PythonInterp 3.5 REQUIRED)
if (PYTHONINTERP_FOUND)
  find_package(PythonLibs "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}" REQUIRED)
endif()

# WITH_COREIR is a macro in the Cython files
# Needs to be set either way
if (WITH_COREIR OR WITH_COREIR_EXTERN)
  set(CYTHON_FLAGS "--compile-time-env WITH_COREIR=ON ${CYTHON_FLAGS}"
      CACHE STRING "Extra flags to the cython compiler." FORCE)
else()
  set(CYTHON_FLAGS "--compile-time-env WITH_COREIR=OFF ${CYTHON_FLAGS}"
      CACHE STRING "Extra flags to the cython compiler." FORCE)
endif()

if (WITH_MSAT_IC3IA)
  set(CYTHON_FLAGS "--compile-time-env WITH_MSAT_IC3IA=ON ${CYTHON_FLAGS}"
      CACHE STRING "Extra flags to the cython compiler." FORCE)
else()
  set(CYTHON_FLAGS "--compile-time-env WITH_MSAT_IC3IA=OFF ${CYTHON_FLAGS}"
      CACHE STRING "Extra flags to the cython compiler." FORCE)
endif()

include_directories("${PROJECT_SOURCE_DIR}/deps/smt-switch/python/smt_switch")

include(FindPythonExtensions)

find_package(Cython 0.29 REQUIRED)

add_cython_target(pono CXX)

add_library(pono MODULE ${pono})

target_link_libraries(pono pono-lib)
target_include_directories(pono PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(pono)

configure_file(setup.py.in setup.py)

install(TARGETS pono DESTINATION lib)
