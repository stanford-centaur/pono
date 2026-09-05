# Sets up generation of options/version.h, which carries the declared release
# version plus the commit provenance of the build. Expects
# PONO_RELEASE_VERSION to be set already.

find_package(Git QUIET)

# Private to pono-lib: the install rules copy *.h out of the source tree, so
# this generated header is deliberately not installed, and
# options/options.cpp is its only includer.
set(PONO_VERSION_HEADER "${CMAKE_CURRENT_BINARY_DIR}/options/version.h")

set(
  PONO_VERSION_GEN_COMMAND
  "${CMAKE_COMMAND}"
  "-DGIT_EXECUTABLE=${GIT_EXECUTABLE}"
  "-DPONO_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
  "-DPONO_RELEASE_VERSION=${PONO_RELEASE_VERSION}"
  "-DINPUT=${PROJECT_SOURCE_DIR}/options/version.h.in"
  "-DOUTPUT=${PONO_VERSION_HEADER}"
  -P
  "${PROJECT_SOURCE_DIR}/cmake/WritePonoVersionHeader.cmake"
)

# Seed the header so that clangd and other compile_commands.json consumers
# find it before the first build.
execute_process(COMMAND ${PONO_VERSION_GEN_COMMAND})

# Re-derived on every build so the provenance cannot go stale between
# configures. This has to be an always-run target rather than a custom command
# with DEPENDS: --dirty flips when a tracked source is edited, which touches
# no file a command could depend on. configure_file rewrites the header only
# when the version actually changed, so the rebuild stays proportionate.
add_custom_target(
  pono-version-header
  BYPRODUCTS "${PONO_VERSION_HEADER}"
  COMMAND ${PONO_VERSION_GEN_COMMAND}
  COMMENT "Checking pono build provenance"
  VERBATIM
)
