#[=======================================================================[.rst:
FindIC3IA
---------

Finds ic3ia (https://es-static.fbk.eu/people/griggio/ic3ia/index.html), an
implementation of IC3 modulo theories with implicit predicate abstraction.

ic3ia has no install step -- ``contrib/setup-ic3ia.sh`` builds it in place --
so this module is the one here that cannot rely on ``IC3IA_ROOT`` alone:
headers stay in the source tree and the library in its build directory,
neither of which is a location ``find_package()`` searches under a prefix.
Both lookups therefore need explicit ``HINTS``.

Hints
^^^^^

``IC3IA_ROOT``
  An ic3ia source tree that has been built, i.e. one containing ``ic3.h`` and
  ``build/libic3ia.a``. ``contrib/setup-ic3ia.sh`` produces one in
  ``deps/ic3ia``.

Imported Targets
^^^^^^^^^^^^^^^^^

``IC3IA::ic3ia``
  The ic3ia library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``IC3IA_FOUND``
  True if the ic3ia library and headers were found.
#]=======================================================================]

# These are cached (as find_path()/find_library() results always are), but a
# cache hit from an earlier configure with a different IC3IA_ROOT would
# otherwise stick around unsearched.
#
# pono includes ic3ia's headers as e.g. "ic3ia/ic3.h" -- it cannot include
# them directly, because some names clash with pono's own (ic3.h). So the
# include directory is the *parent* of the ic3ia tree.
unset(IC3IA_INCLUDE_DIR CACHE)
find_path(IC3IA_INCLUDE_DIR NAMES ic3ia/ic3.h HINTS "${IC3IA_ROOT}/..")

unset(IC3IA_LIBRARY CACHE)
find_library(IC3IA_LIBRARY NAMES ic3ia HINTS "${IC3IA_ROOT}/build")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  IC3IA
  REQUIRED_VARS IC3IA_LIBRARY IC3IA_INCLUDE_DIR
  REASON_FAILURE_MESSAGE
    "Try running contrib/setup-ic3ia.sh, or set --ic3ia-dir to a built ic3ia source tree."
)

if(IC3IA_FOUND AND NOT TARGET IC3IA::ic3ia)
  add_library(IC3IA::ic3ia UNKNOWN IMPORTED)
  set_target_properties(
    IC3IA::ic3ia
    PROPERTIES
      IMPORTED_LOCATION "${IC3IA_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${IC3IA_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(IC3IA_INCLUDE_DIR IC3IA_LIBRARY)
