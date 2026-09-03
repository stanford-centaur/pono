#[=======================================================================[.rst:
FindBtor2Tools
--------------

Finds Btor2Tools (https://github.com/hwmcc/btor2tools), the BTOR2 parser
library. Btor2Tools ships no CMake package config files, so this module
locates its header and library directly.

Hints
^^^^^

``Btor2Tools_ROOT``
  A Btor2Tools install prefix, i.e. the directory passed to its ``--prefix``.
  ``contrib/setup-btor2tools.sh`` installs one into ``deps/install``. Being
  CMake's standard install prefix variable, it needs no explicit ``HINTS``
  below: ``find_package()`` searches it ahead of ``CMAKE_PREFIX_PATH`` and
  the system paths.

Imported Targets
^^^^^^^^^^^^^^^^^

``Btor2Tools::btor2parser``
  The btor2parser library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``Btor2Tools_FOUND``
  True if btor2parser was found.
#]=======================================================================]

# These are cached (as find_path()/find_library() results always are), but a
# cache hit from an earlier configure with a different Btor2Tools_ROOT would
# otherwise stick around unsearched.
unset(Btor2Tools_INCLUDE_DIR CACHE)
find_path(Btor2Tools_INCLUDE_DIR NAMES btor2parser.h)

unset(Btor2Tools_LIBRARY CACHE)
find_library(Btor2Tools_LIBRARY NAMES btor2parser)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Btor2Tools
  REQUIRED_VARS Btor2Tools_LIBRARY Btor2Tools_INCLUDE_DIR
  REASON_FAILURE_MESSAGE
    "Try running contrib/setup-btor2tools.sh, or set --btor2tools-dir to an install prefix."
)

if(Btor2Tools_FOUND AND NOT TARGET Btor2Tools::btor2parser)
  add_library(Btor2Tools::btor2parser UNKNOWN IMPORTED)
  set_target_properties(
    Btor2Tools::btor2parser
    PROPERTIES
      IMPORTED_LOCATION "${Btor2Tools_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${Btor2Tools_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(Btor2Tools_INCLUDE_DIR Btor2Tools_LIBRARY)
