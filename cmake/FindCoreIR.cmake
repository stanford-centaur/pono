#[=======================================================================[.rst:
FindCoreIR
----------

Finds CoreIR (https://github.com/rdaly525/coreir) and the verilogAST library
built alongside it. CoreIR ships no CMake package config files, so this module
locates them directly. Both libraries are required: pono's CoreIR frontend
links against each.

CoreIR is only available as a shared library, so it is incompatible with
``PONO_STATIC_EXEC=YES``.

Hints
^^^^^

``CoreIR_ROOT``
  A CoreIR install prefix. ``contrib/setup-coreir.sh`` installs one into
  ``deps/coreir/local``; a system-wide installation is typically ``/usr`` or
  ``/usr/local``. Being CMake's standard install prefix variable, it needs no
  explicit ``HINTS`` below: ``find_package()`` searches it ahead of
  ``CMAKE_PREFIX_PATH`` and the system paths.

Imported Targets
^^^^^^^^^^^^^^^^^

``CoreIR::coreir``
  The CoreIR library, if found.

``CoreIR::verilogAST``
  The verilogAST library, if found.

Result Variables
^^^^^^^^^^^^^^^^^

``CoreIR_FOUND``
  True if both libraries and the headers were found.
#]=======================================================================]

# These are cached (as find_path()/find_library() results always are), but a
# cache hit from an earlier configure with a different CoreIR_ROOT would
# otherwise stick around unsearched.
unset(CoreIR_INCLUDE_DIR CACHE)
find_path(CoreIR_INCLUDE_DIR NAMES coreir.h)

unset(CoreIR_LIBRARY CACHE)
find_library(CoreIR_LIBRARY NAMES coreir)

unset(CoreIR_verilogAST_LIBRARY CACHE)
find_library(CoreIR_verilogAST_LIBRARY NAMES verilogAST)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CoreIR
  REQUIRED_VARS CoreIR_LIBRARY CoreIR_verilogAST_LIBRARY CoreIR_INCLUDE_DIR
  REASON_FAILURE_MESSAGE
    "Try running contrib/setup-coreir.sh, or set --coreir-dir to an install prefix."
)

if(CoreIR_FOUND)
  if(NOT TARGET CoreIR::coreir)
    add_library(CoreIR::coreir UNKNOWN IMPORTED)
    set_target_properties(
      CoreIR::coreir
      PROPERTIES
        IMPORTED_LOCATION "${CoreIR_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CoreIR_INCLUDE_DIR}"
    )
  endif()
  if(NOT TARGET CoreIR::verilogAST)
    add_library(CoreIR::verilogAST UNKNOWN IMPORTED)
    set_target_properties(
      CoreIR::verilogAST
      PROPERTIES
        IMPORTED_LOCATION "${CoreIR_verilogAST_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${CoreIR_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(CoreIR_INCLUDE_DIR CoreIR_LIBRARY CoreIR_verilogAST_LIBRARY)
