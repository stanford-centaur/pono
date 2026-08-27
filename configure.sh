#!/usr/bin/env bash
# Syntax and structure borrowed from cvc5's configure.sh script
#
# The flags listed in usage(), the variable defaults below it, and the
# case statement in the parsing loop should be kept grouped by type
# (directory settings, build flags, optional features/backends), then
# alphabetically within each group, and in sync with each other.

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

die() {
  echo "*** $0: $*" 1>&2
  exit 1
}

root_dir=$(pwd)

# Directory settings
build_dir=$root_dir/build
install_prefix=/usr/local
smt_switch_dir=$root_dir/deps/smt-switch

# Build flags
build_type=Release
static_exec=NO
lib_type=SHARED

# Optional features / backends
docs=OFF
system_gtest=ON
python=OFF
with_boolector=OFF
with_coreir=OFF
with_coreir_extern=OFF
with_msat=OFF
with_msat_ic3ia=OFF
with_profiling=OFF
with_yices2=OFF
with_z3=OFF

while [[ $# -gt 0 ]]; do
  case $1 in
    -h | --help) usage ;;

    # Directory settings
    --build-dir) die "missing argument to $1 (see -h)" ;;
    --build-dir=*) build_dir=${1##*=} ;;
    --prefix) die "missing argument to $1 (see -h)" ;;
    --prefix=*) install_prefix=${1##*=} ;;
    --smt-switch-dir) die "missing argument to $1 (see -h)" ;;
    --smt-switch-dir=*) smt_switch_dir=${1##*=} ;;

    # Build flags
    --debug)
      build_type=Debug
      ;;
    --static)
      static_exec=YES
      lib_type=STATIC
      ;;
    --static-lib)
      lib_type=STATIC
      ;;

    # Optional features / backends
    --docs) docs=ON ;;
    --no-system-gtest) system_gtest=OFF ;;
    --python)
      python=ON
      ;;
    --with-btor) with_boolector=ON ;;
    --with-coreir) with_coreir=ON ;;
    --with-coreir-extern) with_coreir_extern=ON ;;
    --with-msat) with_msat=ON ;;
    --with-msat-ic3ia) with_msat_ic3ia=ON ;;
    --with-profiling) with_profiling=ON ;;
    --with-yices2) with_yices2=ON ;;
    --with-z3) with_z3=ON ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

[[ $lib_type == STATIC ]] &&
  { [[ $with_coreir == ON ]] || [[ $with_coreir_extern == ON ]]; } &&
  die "CoreIR and static build are incompatible, must omit either" \
    "'--static/--static-lib' or '--with-coreir/--with-coreir-extern'"

# Every option below is passed unconditionally, not just when it differs from
# its default, so that re-running configure.sh on an existing build
# directory always overwrites CMake's cache with the currently requested
# value (e.g. omitting --with-z3 on a later run correctly turns it back off,
# instead of leaving an earlier run's ON stuck in the cache). That in turn
# means we don't have to wipe the build directory to get a correct
# reconfigure.
#
# These are plain CMake variable names, not the bash flags above (which
# don't all map 1:1, e.g. --with-btor -> WITH_BOOLECTOR), so they're ordered
# purely alphabetically rather than matching the flag grouping.
cmake_opts=(
  -DBUILD_DOCS="$docs"
  -DBUILD_PYTHON_BINDINGS="$python"
  -DCMAKE_BUILD_TYPE="$build_type"
  -DCMAKE_INSTALL_PREFIX="$install_prefix"
  -DPONO_LIB_TYPE="$lib_type"
  -DPONO_STATIC_EXEC="$static_exec"
  -DSMT_SWITCH_DIR="$smt_switch_dir"
  -DSYSTEM_GTEST="$system_gtest"
  -DWITH_BOOLECTOR="$with_boolector"
  -DWITH_COREIR="$with_coreir"
  -DWITH_COREIR_EXTERN="$with_coreir_extern"
  -DWITH_MSAT="$with_msat"
  -DWITH_MSAT_IC3IA="$with_msat_ic3ia"
  -DWITH_PROFILING="$with_profiling"
  -DWITH_YICES2="$with_yices2"
  -DWITH_Z3="$with_z3"
)

echo "Running with cmake options: ${cmake_opts[*]}"

cmake -S "$root_dir" -B "$build_dir" "${cmake_opts[@]}" 2>&1
