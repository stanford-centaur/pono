#!/usr/bin/env bash
# The flags listed in usage(), the cmake_vars defaults below it, and the
# case statement in the parsing loop should be kept grouped by type
# (directory settings, build flags, optional features/backends), then
# alphabetically within each group, and in sync with each other.

die() {
  echo "*** $0: $*" 1>&2
  exit 1
}

# Associative arrays (cmake_vars, below) require bash 4+. macOS ships bash
# 3.2 as /bin/bash (Apple won't ship GPLv3 bash), so macOS users need a
# newer bash from Homebrew -- see README.md.
if ((BASH_VERSINFO[0] < 4)); then
  die "requires bash >= 4, but running under bash $BASH_VERSION (on macOS:" \
    "brew install bash, then re-run with that bash)"
fi

usage() {
  cat <<EOF
Usage: $0 [<option> ...]

Configures the CMake build environment.

-h, --help              display this message and exit

Directory settings:
--build-dir=STR         custom build directory (default: build)
--prefix=STR            install directory (default: /usr/local)
--smt-switch-dir=STR    custom smt-switch directory (default: deps/smt-switch)

Build flags:
--debug                 disable optimizations and include debug symbols (default: release build)
--static                build a static executable (default: dynamic); implies --static-lib
--static-lib            build a static library (default: shared)

Optional features / backends (default: off for all):
--docs                  build HTML documentation with Doxygen
--no-system-gtest       do not use system GTest sources; forces download
--python                compile with Python bindings
--with-btor             build with Boolector
--with-coreir           build the CoreIR frontend
--with-coreir-extern    build the CoreIR frontend using an installation in /usr/local/lib
--with-msat             build with MathSAT which has a custom non-BSD compliant license
--with-msat-ic3ia       build with the open-source IC3IA implementation as a backend
--with-profiling        build with gperftools for profiling
--with-yices2           build with Yices2 which has a custom non-BSD compliant license
--with-z3               build with Z3
EOF
  exit 0
}

root_dir=$(pwd)
build_dir=$root_dir/build
declare -A cmake_vars=(
  # Directory settings
  [CMAKE_INSTALL_PREFIX]=/usr/local
  [SMT_SWITCH_DIR]="$root_dir/deps/smt-switch"

  # Build flags
  [CMAKE_BUILD_TYPE]=Release
  [PONO_LIB_TYPE]=SHARED
  [PONO_STATIC_EXEC]=NO

  # Optional features / backends
  [BUILD_DOCS]=OFF
  [BUILD_PYTHON_BINDINGS]=OFF
  [SYSTEM_GTEST]=ON
  [WITH_BOOLECTOR]=OFF
  [WITH_COREIR]=OFF
  [WITH_COREIR_EXTERN]=OFF
  [WITH_MSAT]=OFF
  [WITH_MSAT_IC3IA]=OFF
  [WITH_PROFILING]=OFF
  [WITH_YICES2]=OFF
  [WITH_Z3]=OFF
)

while [[ $# -gt 0 ]]; do
  case $1 in
    -h | --help) usage ;;

    # Directory settings
    --build-dir) die "missing argument to $1 (see -h)" ;;
    --build-dir=*) build_dir=${1##*=} ;;
    --prefix) die "missing argument to $1 (see -h)" ;;
    --prefix=*) cmake_vars[CMAKE_INSTALL_PREFIX]=${1##*=} ;;
    --smt-switch-dir) die "missing argument to $1 (see -h)" ;;
    --smt-switch-dir=*) cmake_vars[SMT_SWITCH_DIR]=${1##*=} ;;

    # Build flags
    --debug) cmake_vars[CMAKE_BUILD_TYPE]=Debug ;;
    # ;& falls through, so --static also runs --static-lib's action below.
    --static) cmake_vars[PONO_STATIC_EXEC]=YES ;&
    --static-lib) cmake_vars[PONO_LIB_TYPE]=STATIC ;;

    # Optional features / backends
    --docs) cmake_vars[BUILD_DOCS]=ON ;;
    --no-system-gtest) cmake_vars[SYSTEM_GTEST]=OFF ;;
    --python) cmake_vars[BUILD_PYTHON_BINDINGS]=ON ;;
    --with-btor) cmake_vars[WITH_BOOLECTOR]=ON ;;
    --with-coreir) cmake_vars[WITH_COREIR]=ON ;;
    --with-coreir-extern) cmake_vars[WITH_COREIR_EXTERN]=ON ;;
    --with-msat) cmake_vars[WITH_MSAT]=ON ;;
    --with-msat-ic3ia) cmake_vars[WITH_MSAT_IC3IA]=ON ;;
    --with-profiling) cmake_vars[WITH_PROFILING]=ON ;;
    --with-yices2) cmake_vars[WITH_YICES2]=ON ;;
    --with-z3) cmake_vars[WITH_Z3]=ON ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

[[ ${cmake_vars[PONO_LIB_TYPE]} == STATIC ]] &&
  { [[ ${cmake_vars[WITH_COREIR]} == ON ]] ||
    [[ ${cmake_vars[WITH_COREIR_EXTERN]} == ON ]]; } &&
  die "CoreIR and static build are incompatible, must omit either" \
    "'--static/--static-lib' or '--with-coreir/--with-coreir-extern'"

# Every option below is passed unconditionally, not just when it differs from
# its default, so that re-running configure.sh on an existing build
# directory always overwrites CMake's cache with the currently requested
# value (e.g. omitting --with-z3 on a later run correctly turns it back off,
# instead of leaving an earlier run's ON stuck in the cache).
cmake_opts=()
for cmake_var in "${!cmake_vars[@]}"; do
  cmake_opts+=("-D$cmake_var=${cmake_vars[$cmake_var]}")
done
cmake -S "$root_dir" -B "$build_dir" "${cmake_opts[@]}" 2>&1
