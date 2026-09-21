# Provisioning for fmt, an implementation detail of utils/logger.h, so there
# is no configure flag pointing at it: use a system copy when one is new
# enough, otherwise fetch one so that building pono needs no extra setup
# step.
#
# Defines fmt::fmt-header-only, the target pono-lib links, and reads
# PONO_LINK_STATIC to decide whether a system copy is usable.

# A fully static link needs that copy to be a static library. Nothing
# links the compiled fmt::fmt yet, since pono uses only the header-only
# target, but a shared-only system copy would break -static for the first
# consumer that does. find_package() is a one-way door here: the imported
# fmt::fmt it defines collides with the ALIAS of that name in fmt's own
# CMakeLists.txt, so a fetch is no longer possible once it has succeeded.
# find_library() defines no targets, so which kind is installed can be
# settled before then.
#
# Cached as every find_library() result is, so drop the entry first: the
# answer depends on PONO_LINK_STATIC as well as on what is installed, and
# either can change between configures of one build directory.
if(PONO_LINK_STATIC)
  unset(PONO_FMT_STATIC_LIBRARY CACHE)
  find_library(PONO_FMT_STATIC_LIBRARY NAMES libfmt.a)
  mark_as_advanced(PONO_FMT_STATIC_LIBRARY)
endif()

if(NOT PONO_LINK_STATIC OR PONO_FMT_STATIC_LIBRARY)
  find_package(fmt 9 QUIET)
endif()

if(fmt_FOUND)
  if(PONO_LINK_STATIC)
    # find_library() and find_package() rank prefixes independently, so the
    # archive found above need not be the one this package describes.
    get_target_property(fmt_library_type fmt::fmt TYPE)
    if(fmt_library_type STREQUAL "SHARED_LIBRARY")
      message(
        FATAL_ERROR
        "The fmt at ${fmt_DIR} provides only a shared library, which "
        "PONO_STATIC_EXEC=YES cannot link. Put the prefix owning "
        "${PONO_FMT_STATIC_LIBRARY} first in CMAKE_PREFIX_PATH, or configure "
        "with -DCMAKE_DISABLE_FIND_PACKAGE_fmt=ON to build a fetched fmt."
      )
    endif()
  endif()
else()
  # The installed utils/logger.h includes <fmt/format.h>, so a prefix built
  # against a fetched fmt has to carry that fmt along to be usable. fmt
  # defaults this off for a subproject build.
  set(FMT_INSTALL ON)
  # fmt takes its library kind from BUILD_SHARED_LIBS, and a shared one would
  # hand a static pono-bin the same unlinkable fmt this fetch exists to
  # avoid. The normal variable shadows any cache entry until it is put back,
  # which has to happen before tests/, where googletest reads the same
  # variable.
  set(pono_build_shared_libs "${BUILD_SHARED_LIBS}")
  set(BUILD_SHARED_LIBS OFF)
  include(FetchContent)
  FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt
    GIT_TAG 1be298e1bd68957e4cd352e1f676f00e07dcfb57 # 12.2.0
  )
  FetchContent_MakeAvailable(fmt)
  set(BUILD_SHARED_LIBS "${pono_build_shared_libs}")
  unset(pono_build_shared_libs)
  # A fetched fmt is built in-tree, so its headers are not treated as system
  # headers the way an imported target's are, and its deprecations would be
  # reported against every file including utils/logger.h. Set the property on
  # the real target since properties cannot be set through an alias.
  get_target_property(fmt_include_dirs fmt-header-only INTERFACE_INCLUDE_DIRECTORIES)
  set_target_properties(
    fmt-header-only
    PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${fmt_include_dirs}"
  )
endif()
