#[=======================================================================[.rst:
FindSmtSwitch
-------------

Finds an already-built smt-switch installation (https://github.com/stanford-centaur/smt-switch),
as produced by its own ``./configure.sh --prefix=... && cmake --build . && cmake --install .``
workflow (this is what ``contrib/setup-smt-switch.sh`` runs).

smt-switch does not ship CMake package config files, so this module locates
its headers/libraries directly and wraps them in imported targets. The
target names deliberately match smt-switch's own (unnamespaced) CMake
target names -- ``smt-switch``, ``smt-switch-bitwuzla``, etc. -- rather than
using a ``SmtSwitch::`` namespace, so that a future in-tree build (e.g. via
``add_subdirectory()`` on a FetchContent-fetched checkout, which defines
targets with these exact names) can be swapped in without changing any
``target_link_libraries()`` call sites.

Hints
^^^^^

``SmtSwitch_ROOT``
  An smt-switch install prefix, i.e. the directory passed to smt-switch's own
  ``--prefix``. Headers and libraries are expected directly under
  ``${SmtSwitch_ROOT}/include`` and ``${SmtSwitch_ROOT}/lib``. Being CMake's
  standard install prefix variable, it needs no explicit ``HINTS`` below:
  ``find_package()`` searches it ahead of ``CMAKE_PREFIX_PATH`` and the
  system paths.

``Bitwuzla_ROOT``
  A bitwuzla install prefix, containing ``lib/pkgconfig/bitwuzla.pc``.
  smt-switch's own ``contrib/setup-bitwuzla.sh`` installs one into
  ``<smt-switch checkout>/deps/install``. Only consulted when the
  ``bitwuzla`` component is requested, and applied by hand rather than by
  ``find_package()``, since bitwuzla is located with pkg-config.

Requesting the ``z3`` component additionally requires a Z3 installation
discoverable via its own CMake package config, since smt-switch's own Z3
backend links against it directly. Set ``Z3_ROOT`` to select one;
smt-switch's own ``contrib/setup-z3.sh`` installs Z3 into
``<smt-switch checkout>/deps/install``. ``Z3_ROOT`` is searched before any
system-wide Z3 installation.

Components
^^^^^^^^^^

Each of ``bitwuzla``, ``cvc5``, ``btor``, ``msat``, ``yices2``, ``z3`` may be
requested as a component. Only requested components are searched for and
given imported targets.

Imported Targets
^^^^^^^^^^^^^^^^

``smt-switch``
  The core smt-switch library, always defined if found.

``smt-switch-<component>``
  One per requested/found component, e.g. ``smt-switch-bitwuzla``.

Result Variables
^^^^^^^^^^^^^^^^^

``SmtSwitch_FOUND``
  True if the core library and all requested components were found.

``SmtSwitch_INCLUDE_DIR``
  Directory containing the ``smt-switch/`` header directory (i.e. the
  directory to add to the include path so that ``#include
  "smt-switch/smt.h"`` resolves).
#]=======================================================================]

# These are cached (as find_path()/find_library() results always are), but
# a cache hit from an earlier configure with a different SmtSwitch_ROOT
# would otherwise stick around unsearched. Unsetting first forces a fresh
# search against the current SmtSwitch_ROOT on every configure.
unset(SmtSwitch_INCLUDE_DIR CACHE)
find_path(SmtSwitch_INCLUDE_DIR NAMES smt-switch/smt.h)

unset(SmtSwitch_LIBRARY CACHE)
find_library(SmtSwitch_LIBRARY NAMES smt-switch)

set(_SmtSwitch_required_vars SmtSwitch_LIBRARY SmtSwitch_INCLUDE_DIR)

# Components that need extra transitive link/include requirements beyond
# just linking the core `smt-switch` target are handled individually below.
foreach(_comp ${SmtSwitch_FIND_COMPONENTS})
  unset(SmtSwitch_${_comp}_LIBRARY CACHE)
  find_library(SmtSwitch_${_comp}_LIBRARY NAMES smt-switch-${_comp})
  if(SmtSwitch_${_comp}_LIBRARY)
    set(SmtSwitch_${_comp}_FOUND TRUE)
  else()
    set(SmtSwitch_${_comp}_FOUND FALSE)
  endif()
  list(APPEND _SmtSwitch_required_vars SmtSwitch_${_comp}_LIBRARY)
endforeach()

# SmtSwitch_ROOT used to be SMT_SWITCH_DIR and took the root of a built
# smt-switch checkout rather than an install prefix, so catch a value that
# still looks like one instead of reporting it as simply not found.
if(NOT SmtSwitch_INCLUDE_DIR AND EXISTS "${SmtSwitch_ROOT}/local/include/smt-switch/smt.h")
  message(
    FATAL_ERROR
    "No smt-switch installation at ${SmtSwitch_ROOT}, but one exists at "
    "${SmtSwitch_ROOT}/local. This now takes an install prefix rather than "
    "the checkout it was built from. Try --smt-switch-dir=${SmtSwitch_ROOT}/local"
  )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  SmtSwitch
  REQUIRED_VARS ${_SmtSwitch_required_vars}
  HANDLE_COMPONENTS
)

if(SmtSwitch_FOUND AND NOT TARGET smt-switch)
  add_library(smt-switch STATIC IMPORTED GLOBAL)
  set_target_properties(
    smt-switch
    PROPERTIES
      IMPORTED_LOCATION "${SmtSwitch_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${SmtSwitch_INCLUDE_DIR}"
  )
  find_package(Threads REQUIRED)
  set_property(TARGET smt-switch APPEND PROPERTY INTERFACE_LINK_LIBRARIES Threads::Threads)

  if("bitwuzla" IN_LIST SmtSwitch_FIND_COMPONENTS AND NOT TARGET smt-switch-bitwuzla)
    add_library(smt-switch-bitwuzla STATIC IMPORTED GLOBAL)
    set_target_properties(
      smt-switch-bitwuzla
      PROPERTIES IMPORTED_LOCATION "${SmtSwitch_bitwuzla_LIBRARY}"
    )
    find_package(PkgConfig REQUIRED)
    # pkg_check_modules() is not a find_*() command, so it does not consult
    # Bitwuzla_ROOT on its own; it does search CMAKE_PREFIX_PATH.
    if(Bitwuzla_ROOT)
      list(APPEND CMAKE_PREFIX_PATH "${Bitwuzla_ROOT}")
    endif()
    # Deliberately not IMPORTED_TARGET: that mode resolves each transitive
    # library name (e.g. mpfr, a public "Requires:" of bitwuzla.pc) via
    # find_library(), which prefers .so over .a on Linux. That breaks fully
    # static executables (PONO_STATIC_EXEC=YES) even when a static archive
    # is available, since the linker is handed an explicit .so path instead
    # of a bare -l flag it could resolve itself. Using the raw pkg-config
    # flags string avoids that.
    pkg_check_modules(SMTSWITCH_BITWUZLA bitwuzla)
    if(NOT SMTSWITCH_BITWUZLA_FOUND)
      message(
        FATAL_ERROR
        "Could not find bitwuzla.pc. Set --bitwuzla-dir to a bitwuzla install "
        "prefix containing lib/pkgconfig/bitwuzla.pc (searched: ${Bitwuzla_ROOT})"
      )
    endif()
    set_property(
      TARGET smt-switch-bitwuzla
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SMTSWITCH_BITWUZLA_INCLUDE_DIRS}
    )
    set_property(
      TARGET smt-switch-bitwuzla
      APPEND
      PROPERTY INTERFACE_LINK_LIBRARIES smt-switch "${SMTSWITCH_BITWUZLA_LDFLAGS}"
    )
  endif()

  if("z3" IN_LIST SmtSwitch_FIND_COMPONENTS AND NOT TARGET smt-switch-z3)
    add_library(smt-switch-z3 STATIC IMPORTED GLOBAL)
    set_target_properties(smt-switch-z3 PROPERTIES IMPORTED_LOCATION "${SmtSwitch_z3_LIBRARY}")
    # Z3_ROOT is searched ahead of CMAKE_PREFIX_PATH and the system paths.
    # Clear the cached config directory first, or an earlier configure's Z3
    # would survive a change of Z3_ROOT.
    unset(Z3_DIR CACHE)
    find_package(Z3 QUIET)
    if(NOT Z3_FOUND)
      message(
        FATAL_ERROR
        "Could not find Z3Config.cmake. Set --z3-dir to a Z3 install prefix "
        "containing lib/cmake/z3 (searched: ${Z3_ROOT})"
      )
    endif()
    set_property(
      TARGET smt-switch-z3
      APPEND
      PROPERTY INTERFACE_LINK_LIBRARIES smt-switch z3::libz3
    )
  endif()

  foreach(_comp cvc5 btor msat yices2)
    if(_comp IN_LIST SmtSwitch_FIND_COMPONENTS AND NOT TARGET smt-switch-${_comp})
      add_library(smt-switch-${_comp} STATIC IMPORTED GLOBAL)
      set_target_properties(
        smt-switch-${_comp}
        PROPERTIES IMPORTED_LOCATION "${SmtSwitch_${_comp}_LIBRARY}"
      )
      set_property(TARGET smt-switch-${_comp} APPEND PROPERTY INTERFACE_LINK_LIBRARIES smt-switch)
    endif()
  endforeach()
endif()

mark_as_advanced(SmtSwitch_INCLUDE_DIR SmtSwitch_LIBRARY)
