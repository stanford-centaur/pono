#[=======================================================================[.rst:
FindSmtSwitch
-------------

Finds an already-built smt-switch installation (https://github.com/stanford-centaur/smt-switch),
as produced by its own ``./configure.sh --prefix=local ... && cmake --build . && cmake --install .``
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

``SMT_SWITCH_DIR``
  Root of a built smt-switch checkout, i.e. the directory passed as
  ``--prefix=local`` was run from. Headers/libraries are expected under
  ``${SMT_SWITCH_DIR}/local``. Bitwuzla's pkg-config files (built by
  smt-switch's own ``contrib/setup-bitwuzla.sh``) are expected under
  ``${SMT_SWITCH_DIR}/deps/install``.

Requesting the ``z3`` component additionally requires a Z3 installation
discoverable via its own CMake package config, since smt-switch's own Z3
backend links against it directly. smt-switch's own ``contrib/setup-z3.sh``
installs Z3 into ``${SMT_SWITCH_DIR}/deps/install``; that location is always
preferred over a system-wide Z3 installation if both are present.

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

if(NOT DEFINED SMT_SWITCH_DIR)
  set(SMT_SWITCH_DIR "${PROJECT_SOURCE_DIR}/deps/smt-switch")
endif()

find_path(SmtSwitch_INCLUDE_DIR NAMES smt-switch/smt.h HINTS "${SMT_SWITCH_DIR}/local/include")

find_library(SmtSwitch_LIBRARY NAMES smt-switch HINTS "${SMT_SWITCH_DIR}/local/lib")

set(_SmtSwitch_required_vars SmtSwitch_LIBRARY SmtSwitch_INCLUDE_DIR)

# Components that need extra transitive link/include requirements beyond
# just linking the core `smt-switch` target are handled individually below.
foreach(_comp ${SmtSwitch_FIND_COMPONENTS})
  find_library(
    SmtSwitch_${_comp}_LIBRARY
    NAMES smt-switch-${_comp}
    HINTS "${SMT_SWITCH_DIR}/local/lib"
  )
  if(SmtSwitch_${_comp}_LIBRARY)
    set(SmtSwitch_${_comp}_FOUND TRUE)
  else()
    set(SmtSwitch_${_comp}_FOUND FALSE)
  endif()
  list(APPEND _SmtSwitch_required_vars SmtSwitch_${_comp}_LIBRARY)
endforeach()

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
    list(APPEND CMAKE_PREFIX_PATH "${SMT_SWITCH_DIR}/deps/install")
    # Deliberately not IMPORTED_TARGET: that mode resolves each transitive
    # library name (e.g. mpfr, a public "Requires:" of bitwuzla.pc) via
    # find_library(), which prefers .so over .a on Linux. That breaks fully
    # static executables (PONO_STATIC_EXEC=YES) even when a static archive
    # is available, since the linker is handed an explicit .so path instead
    # of a bare -l flag it could resolve itself. Using the raw pkg-config
    # flags string avoids that.
    pkg_check_modules(SMTSWITCH_BITWUZLA REQUIRED bitwuzla)
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
    # Prefer the Z3 build smt-switch's own contrib/setup-z3.sh installs into
    # deps/install over any system-wide Z3 installation.
    find_package(Z3 QUIET PATHS "${SMT_SWITCH_DIR}/deps/install" NO_DEFAULT_PATH)
    if(NOT Z3_FOUND)
      find_package(Z3 REQUIRED)
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
