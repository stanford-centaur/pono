#[=======================================================================[.rst:
FindGperftools
--------------

Finds the CPU profiler from gperftools (https://github.com/gperftools/gperftools)
via pkg-config. gperftools ships no CMake package config files, so this module
locates it directly.

Imported Targets
^^^^^^^^^^^^^^^^^

``Gperftools::profiler``
  The libprofiler library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``Gperftools_FOUND``
  True if libprofiler was found.
#]=======================================================================]

find_package(PkgConfig REQUIRED)
pkg_check_modules(Gperftools QUIET libprofiler)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Gperftools
  REQUIRED_VARS Gperftools_LDFLAGS
  VERSION_VAR Gperftools_VERSION
  REASON_FAILURE_MESSAGE
    "No libprofiler.pc found. Install gperftools to build with --with-profiling."
)

if(Gperftools_FOUND AND NOT TARGET Gperftools::profiler)
  add_library(Gperftools::profiler INTERFACE IMPORTED)
  set_target_properties(
    Gperftools::profiler
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Gperftools_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${Gperftools_LDFLAGS}"
  )
endif()
