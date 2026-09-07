#[=======================================================================[.rst:
FindBitwuzla
------------

Finds the Bitwuzla SMT solver (https://github.com/bitwuzla/bitwuzla) via
pkg-config. Bitwuzla is built with Meson and ships no CMake package config
files, so this module locates it directly.

Hints
^^^^^

``Bitwuzla_ROOT``
  A bitwuzla install prefix whose pkgconfig directory holds ``bitwuzla.pc``.
  smt-switch's own ``contrib/setup-bitwuzla.sh`` installs one into
  ``<smt-switch checkout>/deps/install``. If unset, only the system
  pkg-config search path is used.

Imported Targets
^^^^^^^^^^^^^^^^^

``Bitwuzla::bitwuzla``
  The bitwuzla library and its transitive dependencies, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``Bitwuzla_FOUND``
  True if bitwuzla was found.
#]=======================================================================]

find_package(PkgConfig REQUIRED)

# pkg_check_modules() is not a find_*() command, so it ignores Bitwuzla_ROOT.
# It does consult CMAKE_PREFIX_PATH, and works out which pkgconfig
# subdirectory a prefix uses on this platform -- lib/<arch>/pkgconfig on
# Debian derivatives, lib64 or libdata elsewhere, then lib and share. Setting
# PKG_CONFIG_PATH by hand instead would have to hardcode one of those.
#
# Restore the variable afterwards: a find module runs in its caller's
# directory scope, so leaving it extended would let this prefix satisfy every
# later find_*() call there, e.g. supplying a Z3 that overrides Z3_ROOT.
set(_Bitwuzla_saved_prefix_path "${CMAKE_PREFIX_PATH}")
if(Bitwuzla_ROOT)
  list(PREPEND CMAKE_PREFIX_PATH "${Bitwuzla_ROOT}")
endif()
pkg_check_modules(Bitwuzla QUIET bitwuzla)
set(CMAKE_PREFIX_PATH "${_Bitwuzla_saved_prefix_path}")
unset(_Bitwuzla_saved_prefix_path)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Bitwuzla
  REQUIRED_VARS Bitwuzla_INCLUDE_DIRS Bitwuzla_LDFLAGS
  VERSION_VAR Bitwuzla_VERSION
  REASON_FAILURE_MESSAGE "No bitwuzla.pc found. Set --bitwuzla-dir to a bitwuzla install prefix."
)

if(Bitwuzla_FOUND AND NOT TARGET Bitwuzla::bitwuzla)
  add_library(Bitwuzla::bitwuzla INTERFACE IMPORTED)
  # Deliberately not IMPORTED_TARGET on pkg_check_modules(): that mode
  # resolves each transitive library name (e.g. mpfr, a public "Requires:" of
  # bitwuzla.pc) via find_library(), which prefers .so over .a on Linux. That
  # breaks fully static executables (PONO_STATIC_EXEC=YES) even when a static
  # archive is available, since the linker is handed an explicit .so path
  # instead of a bare -l flag it could resolve itself. Using the raw
  # pkg-config flags string avoids that.
  set_target_properties(
    Bitwuzla::bitwuzla
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Bitwuzla_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${Bitwuzla_LDFLAGS}"
  )
endif()
