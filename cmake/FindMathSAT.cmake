#[=======================================================================[.rst:
FindMathSAT
-----------

Finds the MathSAT 5 SMT solver's headers (https://mathsat.fbk.eu). MathSAT is
distributed as a prebuilt archive and ships no CMake package config files, so
this module locates it directly.

Only the headers are located, not ``libmathsat``. pono needs direct access to
MathSAT's C API to set solver options, but links MathSAT through smt-switch:
a static smt-switch build repacks MathSAT's objects into
``libsmt-switch-msat.a``, so also linking the standalone archive would risk
duplicate definitions.

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
  Interface target carrying MathSAT's include directory, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``MathSAT_FOUND``
  True if the MathSAT headers were found.
#]=======================================================================]

# Cached (as find_path() results always are), but a cache hit from an earlier
# configure with a different MathSAT_ROOT would otherwise stick around
# unsearched.
unset(MathSAT_INCLUDE_DIR CACHE)
find_path(MathSAT_INCLUDE_DIR NAMES mathsat.h)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  MathSAT
  REQUIRED_VARS MathSAT_INCLUDE_DIR
  REASON_FAILURE_MESSAGE
    "Try running ci-scripts/setup-msat.sh, or set --msat-dir to an install prefix."
)

if(MathSAT_FOUND AND NOT TARGET MathSAT::mathsat)
  add_library(MathSAT::mathsat INTERFACE IMPORTED)
  set_target_properties(
    MathSAT::mathsat
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MathSAT_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(MathSAT_INCLUDE_DIR)
