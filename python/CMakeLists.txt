
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

# Everything from here up to line 'include(FindPythonExtensions)' is copied
# from smt-switch python/CMakefile:
#   https://github.com/stanford-centaur/smt-switch/blob/499621b009ca0e86a65339aed59157bfaa874776/python/CMakeLists.txt
# Note: using same pattern used in cvc5 for scikit-build / python module checking
# Check if given Python module is installed and raises a FATAL_ERROR error
# if the module cannot be found.
function(check_python_module module)
  execute_process(
    COMMAND
    ${PYTHON_EXECUTABLE} -c "import ${module}"
    RESULT_VARIABLE
      RET_MODULE_TEST
    ERROR_QUIET
  )
  set(module_name ${ARGN})
  if(NOT module_name)
    set(module_name ${module})
  endif()

  if(RET_MODULE_TEST)
    message(FATAL_ERROR
        "Could not find module ${module_name} for Python "
        "version ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}. "
        "Make sure to install ${module_name} for this Python version "
        "via \n`${PYTHON_EXECUTABLE} -m pip install ${module_name}'.\n"
        "Note: You need to have pip installed for this Python version.")
  endif()
endfunction()


# Check that scikit-build is installed
# Provides CMake files for Python bindings
# and used to find correct Python includes
# and libraries
check_python_module("skbuild" "scikit-build")

# Find cmake modules distributed with scikit-build
# They are distributed under <scikit-build directory>/resources/cmake
execute_process(
  COMMAND
    ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; \
      import os; import skbuild; \
      cmake_module_path=os.path.join(os.path.dirname(skbuild.__file__), \
        'resources', 'cmake'); print(cmake_module_path, end='')"
  OUTPUT_VARIABLE
    SKBUILD_CMAKE_MODULE_PATH
  RESULT_VARIABLE
    RET_SKBUILD_CMAKE_MODULE_PATH
)

if (NOT EXISTS ${SKBUILD_CMAKE_MODULE_PATH})
  message(FATAL_ERROR "Expected CMake module path from
                       scikit-build at ${SKBUILD_CMAKE_MODULE_PATH}")
endif()

# Add scikit-build cmake files to cmake module path
# Required for Cython target below
list(APPEND CMAKE_MODULE_PATH ${SKBUILD_CMAKE_MODULE_PATH})

include(FindPythonExtensions)

find_package(Cython 0.29 REQUIRED)

add_cython_target(pono CXX)

add_library(pono MODULE ${pono})

target_link_libraries(pono pono-lib)
target_include_directories(pono PUBLIC "${PROJECT_SOURCE_DIR}/deps/smt-switch/local/include/smt-switch")

python_extension_module(pono)

configure_file(setup.py.in setup.py)

install(TARGETS pono DESTINATION lib)
