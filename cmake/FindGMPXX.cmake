#[=======================================================================[.rst:
FindGMPXX
---------

Finds the GNU Multiple Precision Arithmetic Library's C++ bindings
(gmpxx) via pkg-config. GMP ships no CMake package config files, so this
module locates it directly.

Imported Targets
^^^^^^^^^^^^^^^^^

``GMPXX::gmpxx``
  The gmpxx library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``GMPXX_FOUND``
  True if gmpxx was found.
#]=======================================================================]

find_package(PkgConfig REQUIRED)
pkg_check_modules(GMPXX QUIET gmpxx)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  GMPXX
  REQUIRED_VARS GMPXX_INCLUDE_DIRS GMPXX_LDFLAGS
  VERSION_VAR GMPXX_VERSION
)

if(GMPXX_FOUND AND NOT TARGET GMPXX::gmpxx)
  add_library(GMPXX::gmpxx INTERFACE IMPORTED)
  set_target_properties(
    GMPXX::gmpxx
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${GMPXX_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${GMPXX_LDFLAGS}"
  )
endif()
