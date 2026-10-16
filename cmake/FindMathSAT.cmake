#[=======================================================================[.rst:
FindMathSAT
-----------

Finds the MathSAT 5 SMT solver (https://mathsat.fbk.eu). MathSAT is
distributed as a prebuilt archive and ships no CMake package config files, so
this module locates it directly. pono needs direct access to MathSAT's C API
to set solver options.

smt-switch's package locates MathSAT for its own backend with a module of the
same name, which leaves a ``MathSAT::mathsat`` that already exists alone. The
target defined here is therefore the one both link, and has to carry the
library as well as the headers.

Hints
^^^^^

``MathSAT_ROOT``
  A MathSAT install prefix -- the unpacked release archive, which contains
  ``include`` and ``lib`` directories. ``ci-scripts/setup-msat.sh`` unpacks
  one into ``deps/mathsat``. Being CMake's standard install prefix variable,
  it needs no explicit ``HINTS`` below: ``find_package()`` searches it ahead
  of ``CMAKE_PREFIX_PATH`` and the system paths.

Imported Targets
^^^^^^^^^^^^^^^^^

``MathSAT::mathsat``
  The MathSAT library and its include directory, if found. MathSAT needs GMP,
  which smt-switch's MathSAT backend links.

Result Variables
^^^^^^^^^^^^^^^^^

``MathSAT_FOUND``
  True if the MathSAT headers and library were found.
#]=======================================================================]

# Cached (as find_path() and find_library() results always are), but a cache
# hit from an earlier configure with a different MathSAT_ROOT would otherwise
# stick around unsearched.
unset(MathSAT_INCLUDE_DIR CACHE)
find_path(MathSAT_INCLUDE_DIR NAMES mathsat.h)
unset(MathSAT_LIBRARY CACHE)
find_library(MathSAT_LIBRARY NAMES mathsat)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MathSAT
  REQUIRED_VARS MathSAT_LIBRARY MathSAT_INCLUDE_DIR
  REASON_FAILURE_MESSAGE
    "Try running ci-scripts/setup-msat.sh, or set --msat-dir to an install prefix."
)

if(MathSAT_FOUND AND NOT TARGET MathSAT::mathsat)
  add_library(MathSAT::mathsat UNKNOWN IMPORTED)
  set_target_properties(
    MathSAT::mathsat
    PROPERTIES
      IMPORTED_LOCATION "${MathSAT_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MathSAT_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(MathSAT_INCLUDE_DIR MathSAT_LIBRARY)
