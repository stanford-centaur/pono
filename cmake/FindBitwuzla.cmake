#[=======================================================================[.rst:
FindBitwuzla
------------

Finds the Bitwuzla SMT solver (https://github.com/bitwuzla/bitwuzla) via
pkg-config. Bitwuzla is built with Meson and ships no CMake package config
files, so this module locates it directly.

Hints
^^^^^

``Bitwuzla_ROOT``
  A bitwuzla install prefix, containing ``lib/pkgconfig/bitwuzla.pc``.
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
# Extending CMAKE_PREFIX_PATH would reach it, but a find module runs in its
# caller's directory scope, so that would leak into every later find_*() call
# there and let this prefix satisfy unrelated lookups. Point pkg-config at the
# prefix directly instead, and restore the environment afterwards.
set(_Bitwuzla_saved_pkg_config_path "$ENV{PKG_CONFIG_PATH}")
if(Bitwuzla_ROOT)
  set(ENV{PKG_CONFIG_PATH} "${Bitwuzla_ROOT}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif()
pkg_check_modules(Bitwuzla QUIET bitwuzla)
set(ENV{PKG_CONFIG_PATH} "${_Bitwuzla_saved_pkg_config_path}")
unset(_Bitwuzla_saved_pkg_config_path)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Bitwuzla
  REQUIRED_VARS Bitwuzla_INCLUDE_DIRS Bitwuzla_LDFLAGS
  VERSION_VAR Bitwuzla_VERSION
  REASON_FAILURE_MESSAGE
    "No bitwuzla.pc found. Set --bitwuzla-dir to a prefix containing lib/pkgconfig/bitwuzla.pc."
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
