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

# Need to make sure libraries match the interpreter
# Our CMake files use the latest Python finding cmake modules (as of 2020)
# https://cmake.org/cmake/help/v3.12/module/FindPython.html (specifically the versioned one
#  FindPython3
find_package (Python3 COMPONENTS Interpreter Development)


# However, the Python Extensions from scikit-build still use an old version
# We need to set variables so that FindPythonInterp is not run in FindPythonExtensions
# below, but it looks in the right place for the PythonLibs
set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
set(PYTHONINTERP_FOUND ${Python3_Interpreter_FOUND})
set(PYTHON_VERSION_STRING ${Python3_VERSION})
set(PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
set(PYTHON_VERSION_MINOR ${Python3_VERSION_MINOR})
set(PYTHON_VERSION_PATCH ${Python3_VERSION_PATCH})

if (Python3_Development_FOUND)
  set(PYTHON_LIBRARY "${Python3_LIBRARIES}" CACHE STRING "Python libraries")
  set(PYTHON_INCLUDE_DIR "${Python3_INCLUDE_DIRS}" CACHE STRING "Python include directories")
else()
  message(STATUS "Didn't find python library and include directory -- falling back")
  execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "from distutils import sysconfig as s;
print(s.get_python_inc(plat_specific=True));
print(s.get_config_var('LDVERSION') or s.get_config_var('VERSION'));
"
     RESULT_VARIABLE PY_SUCCESS
     OUTPUT_VARIABLE PY_VALS
     ERROR_VARIABLE  PY_ERR)

  if (NOT PY_SUCCESS MATCHES 0)
    message(FATAL "Cannot locate python library...")
  endif()

  # create a list
  string(REGEX REPLACE ";" "\\\\;" PY_VALS ${PY_VALS})
  string(REGEX REPLACE "\n" ";" PY_VALS ${PY_VALS})

  list(GET PY_VALS 0 PY_INCLUDE_DIR)
  list(GET PY_VALS 1 PY_LIB_SUFFIX)

  set(PYTHON_INCLUDE_DIR "${PY_INCLUDE_DIR}" CACHE STRING "Python include directory")

  # since cmake could not find the library, just set the name/version
  # and let the linker figure it out
  set(PYTHON_LIBRARY python${PY_LIB_SUFFIX} CACHE STRING "Python library")

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

include_directories("${PROJECT_SOURCE_DIR}/deps/smt-switch/python")

include(FindPythonExtensions)

find_package(Cython 0.29 REQUIRED)

add_cython_target(pono CXX)

add_library(pono MODULE ${pono})

target_link_libraries(pono pono-lib)
target_include_directories(pono PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(pono)

configure_file(setup.py.in setup.py)

install(TARGETS pono DESTINATION lib)
