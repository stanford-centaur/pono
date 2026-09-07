# Generates the pono version header. Run via `cmake -P`, both once at
# configure time to seed the header and again on every build, so that the
# commit provenance below is re-derived instead of snapshotted.
#
# Inputs, all passed with -D:
#
# GIT_EXECUTABLE       git binary, or empty when git is unavailable
# PONO_SOURCE_DIR      pono source tree to describe
# PONO_RELEASE_VERSION declared version, reported with or without git
# INPUT                header template to configure
# OUTPUT               header to generate

set(PONO_GIT_DESCRIBE "")

if(GIT_EXECUTABLE)
  # --tags is load-bearing: release tags here are lightweight, and plain
  # `git describe` considers only annotated tags. --match confines the search
  # to release tags, so a competition snapshot tag (hwmcc*) landing on
  # mainline cannot win on distance. -C resolves the gitfile pointer that a
  # linked worktree has in place of a .git directory.
  execute_process(
    COMMAND
      "${GIT_EXECUTABLE}" -C "${PONO_SOURCE_DIR}" describe --tags --always --dirty --match
      "v[0-9]*"
    RESULT_VARIABLE GIT_DESCRIBE_STATUS
    OUTPUT_VARIABLE PONO_GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  # Anything short of success (no git history, as in a tarball or a Docker
  # build context) leaves the version without a provenance suffix.
  if(NOT GIT_DESCRIBE_STATUS EQUAL 0)
    set(PONO_GIT_DESCRIBE "")
  endif()
endif()

if(PONO_GIT_DESCRIBE)
  set(PONO_VERSION "${PONO_RELEASE_VERSION} (${PONO_GIT_DESCRIBE})")
else()
  set(PONO_VERSION "${PONO_RELEASE_VERSION}")
endif()

configure_file("${INPUT}" "${OUTPUT}" @ONLY)
