#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
build_dir=$root_dir/build

# These variables, the flags in usage(), and the case statement below should
# stay in sync, grouped by type (build/install directories, dependency
# locations, build flags, optional features), and alphabetical within each
# group.
#
# There is one cm_-prefixed variable per CMake variable. These should ideally
# be in a single associative array, but macOS's default bash (3.2) doesn't
# support those, so each variable is declared individually instead; the shared
# prefix lets the final -D flag loop enumerate them all via ${!cm_@} without a
# separately-maintained list of names.

# Build and install directories
cm_CMAKE_INSTALL_PREFIX=/usr/local

# Dependency locations, each named <Package>_ROOT after CMake's standard
# install prefix input, using the package name the dependency would most
# likely get as a CMake module. Z3Config.cmake and cmake/FindSmtSwitch.cmake
# let find_package() honor Z3_ROOT and SmtSwitch_ROOT automatically; bitwuzla
# (pkg-config) and MathSAT (a bare include path) have no module yet, so their
# roots are applied by hand until one exists.
cm_Bitwuzla_ROOT="$root_dir/deps/smt-switch/deps/install"
cm_MathSAT_ROOT="$root_dir/deps/mathsat"
cm_SmtSwitch_ROOT="$root_dir/deps/smt-switch/local"
cm_Z3_ROOT="$root_dir/deps/smt-switch/deps/install"

# Build flags
cm_CMAKE_BUILD_TYPE=Release
cm_PONO_LIB_TYPE=SHARED
cm_PONO_STATIC_EXEC=NO

# Optional features / backends
cm_BUILD_DOCS=OFF
cm_BUILD_PYTHON_BINDINGS=OFF
cm_SYSTEM_GTEST=ON
cm_WITH_BOOLECTOR=OFF
cm_WITH_COREIR=OFF
cm_WITH_COREIR_EXTERN=OFF
cm_WITH_MSAT=OFF
cm_WITH_MSAT_IC3IA=OFF
cm_WITH_PROFILING=OFF
cm_WITH_YICES2=OFF
cm_WITH_Z3=OFF

# Utility functions
die() { echo "*** $0: $*" 1>&2 && exit 1; }
lowercase() { printf '%s\n' "$1" | tr '[:upper:]' '[:lower:]'; }

usage() {
  cat <<EOF
Usage: $0 [<option> ...]

Configures the CMake build environment.

-h, --help              display this message and exit

Build and install directories:
--build-dir=STR         custom build directory (default: $build_dir)
--prefix=STR            install directory (default: $cm_CMAKE_INSTALL_PREFIX)

Dependency locations:
--bitwuzla-dir=STR      Bitwuzla install prefix (default: $cm_Bitwuzla_ROOT)
--msat-dir=STR          MathSAT install prefix (default: $cm_MathSAT_ROOT)
--smt-switch-dir=STR    smt-switch install prefix (default: $cm_SmtSwitch_ROOT)
--z3-dir=STR            Z3 install prefix (default: $cm_Z3_ROOT)

Build flags:
--debug                 disable optimizations and include debug symbols (default: $(lowercase "$cm_CMAKE_BUILD_TYPE") build)
--static                build a static executable (default: $(lowercase "$cm_PONO_STATIC_EXEC")); implies --static-lib
--static-lib            build a static library (default: $(lowercase "$cm_PONO_LIB_TYPE"))

Optional features / backends (default: false/off for all):
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

while [[ $# -gt 0 ]]; do
  # shellcheck disable=SC2034 # shellcheck can't trace ${!name} indirection
  case $1 in
    -h | --help) usage ;;

    # Build and install directories
    --build-dir) die "missing argument to $1 (see -h)" ;;
    --build-dir=*) build_dir=${1#*=} ;;
    --prefix) die "missing argument to $1 (see -h)" ;;
    --prefix=*) cm_CMAKE_INSTALL_PREFIX=${1#*=} ;;

    # Dependency locations
    --bitwuzla-dir) die "missing argument to $1 (see -h)" ;;
    --bitwuzla-dir=*) cm_Bitwuzla_ROOT=${1#*=} ;;
    --msat-dir) die "missing argument to $1 (see -h)" ;;
    --msat-dir=*) cm_MathSAT_ROOT=${1#*=} ;;
    --smt-switch-dir) die "missing argument to $1 (see -h)" ;;
    --smt-switch-dir=*) cm_SmtSwitch_ROOT=${1#*=} ;;
    --z3-dir) die "missing argument to $1 (see -h)" ;;
    --z3-dir=*) cm_Z3_ROOT=${1#*=} ;;

    # Build flags
    --debug) cm_CMAKE_BUILD_TYPE=Debug ;;
    --static) cm_PONO_STATIC_EXEC=YES ;;
    --static-lib) cm_PONO_LIB_TYPE=STATIC ;;

    # Optional features / backends
    --docs) cm_BUILD_DOCS=ON ;;
    --no-system-gtest) cm_SYSTEM_GTEST=OFF ;;
    --python) cm_BUILD_PYTHON_BINDINGS=ON ;;
    --with-btor) cm_WITH_BOOLECTOR=ON ;;
    --with-coreir) cm_WITH_COREIR=ON ;;
    --with-coreir-extern) cm_WITH_COREIR_EXTERN=ON ;;
    --with-msat) cm_WITH_MSAT=ON ;;
    --with-msat-ic3ia) cm_WITH_MSAT_IC3IA=ON ;;
    --with-profiling) cm_WITH_PROFILING=ON ;;
    --with-yices2) cm_WITH_YICES2=ON ;;
    --with-z3) cm_WITH_Z3=ON ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

# Every option below is passed unconditionally, not just when it differs from
# its default, so that re-running configure.sh on an existing build
# directory always overwrites CMake's cache with the currently requested
# value (e.g. omitting --with-z3 on a later run correctly turns it back off,
# instead of leaving an earlier run's ON stuck in the cache).
cmake_opts=()
for name in "${!cm_@}"; do
  cmake_opts+=("-D${name#cm_}=${!name}")
done
cmake -S "$root_dir" -B "$build_dir" "${cmake_opts[@]}" 2>&1
