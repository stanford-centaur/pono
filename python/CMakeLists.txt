# Check that Python executable and development headers are available.
find_package(Python 3.5 REQUIRED COMPONENTS Interpreter Development.Module)

# Check that cython >3.0 is available.
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c "import cython; print(cython.__version__)"
  RESULT_VARIABLE cython_not_found
  OUTPUT_VARIABLE cython_version
  ERROR_QUIET
)
if(cython_not_found)
  message(FATAL_ERROR "Could not find Python module cython, install via ${Python_EXECUTABLE} -m pip install Cython.")
endif()
string(REPLACE "." ";" cython_version_list ${cython_version})
list(GET cython_version_list 0 cython_version_major)
if(cython_version_major LESS 3)
  message(FATAL_ERROR "Cython version 3.0.0 or later required; currently installed version is ${cython_version}")
endif()

# Set up Python wheel build.
configure_file(pyproject.toml.in pyproject.toml @ONLY)

# Build the actual extension module.
add_subdirectory(pono)
