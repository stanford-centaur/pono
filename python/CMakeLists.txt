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

# WITH_COREIR is a macro in the Cython files
# Needs to be set either way
if (WITH_COREIR)
  set(CYTHON_FLAGS "--compile-time-env WITH_COREIR=ON ${CYTHON_FLAGS}"
      CACHE STRING "Extra flags to the cython compiler." FORCE)
else()
  set(CYTHON_FLAGS "--compile-time-env WITH_COREIR=OFF ${CYTHON_FLAGS}"
                   CACHE STRING "Extra flags to the cython compiler." FORCE)
endif()

find_package(Cython 0.29 REQUIRED)

include_directories("${PROJECT_SOURCE_DIR}/deps/smt-switch/python")

# Need to make sure libraries match the interpreter
# Our CMake files use the latest Python finding cmake modules (as of 2020)
# https://cmake.org/cmake/help/v3.12/module/FindPython.html (specifically the versioned ones
#  FindPython2 and FindPython3

# However, the Python Extensions from scikit-build still use an old version
# We need to set variables so that FindPythonInterp is not run in FindPythonExtensions
# below, but it looks in the right place for the PythonLibs
if (USE_PYTHON2)
  set(PYTHON_EXECUTABLE ${Python2_EXECUTABLE})
  set(PYTHONINTERP_FOUND ${Python2_Interpreter_FOUND})
  set(PYTHON_VERSION_STRING ${Python2_VERSION})
  set(PYTHON_VERSION_MAJOR ${Python2_VERSION_MAJOR})
  set(PYTHON_VERSION_MINOR ${Python2_VERSION_MINOR})
  set(PYTHON_VERSION_PATCH ${Python2_VERSION_PATCH})
else()
  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
  set(PYTHONINTERP_FOUND ${Python3_Interpreter_FOUND})
  set(PYTHON_VERSION_STRING ${Python3_VERSION})
  set(PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
  set(PYTHON_VERSION_MINOR ${Python3_VERSION_MINOR})
  set(PYTHON_VERSION_PATCH ${Python3_VERSION_PATCH})
endif()

include(FindPythonExtensions)

add_cython_target(pono CXX)

add_library(pono MODULE ${pono})

target_link_libraries(pono pono-lib)
target_include_directories(pono PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(pono)

configure_file(setup.py.in setup.py)

install(TARGETS pono DESTINATION lib)
