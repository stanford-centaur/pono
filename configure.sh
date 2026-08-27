#!/usr/bin/env bash

# Syntax and structure borrowed from cvc5's configure.sh script

usage() {
  cat <<EOF
Usage: $0 [<option> ...]

Configures the CMAKE build environment.

-h, --help              display this message and exit
--prefix=STR            install directory       (default: /usr/local/)
--build-dir=STR         custom build directory  (default: build)
--smt-switch-dir=STR    custom smt-switch directory (default: deps/smt-switch)
--with-btor             build with Boolector  (default: off)
--with-msat             build with MathSAT which has a custom non-BSD compliant license.  (default : off)
                        Required for interpolant based model checking
--with-yices2           build with Yices2 which has a custom non-BSD compliant license (default : off)
--with-z3               build with Z3 (default: off)
--with-msat-ic3ia       build with the open-source IC3IA implementation as a backend. (default: off)
--with-coreir           build the CoreIR frontend (default: off)
--with-coreir-extern    build the CoreIR frontend using an installation of coreir in /usr/local/lib (default: off)
--debug                 build debug with debug symbols (default: off)
--docs                  build HTML documentation with Doxygen (default: off)
--python                compile with python bindings (default: off)
--static-lib            build a static library (default: shared)
--static                build a static executable (default: dynamic); implies --static-lib
--with-profiling        build with gperftools for profiling (default: off)
--no-system-gtest       do not use system GTest sources; forces download (default: off)
EOF
  exit 0
}

die() {
  echo "*** $0: $*" 1>&2
  exit 1
}

root_dir=$(pwd)

build_dir=$root_dir/build
smt_switch_dir=$root_dir/deps/smt-switch
install_prefix=/usr/local
build_type=Release
with_boolector=OFF
with_msat=OFF
with_yices2=OFF
with_z3=OFF
with_msat_ic3ia=OFF
with_coreir=OFF
with_coreir_extern=OFF
docs=OFF
python=OFF
lib_type=SHARED
static_exec=NO
with_profiling=OFF
system_gtest=ON

while [[ $# -gt 0 ]]; do
  case $1 in
    -h | --help) usage ;;
    --prefix) die "missing argument to $1 (see -h)" ;;
    --prefix=*)
      install_prefix=${1##*=}
      # Check if install_prefix is an absolute path and if not, make it
      # absolute.
      case $install_prefix in
        /*) ;;                                         # absolute path
        *) install_prefix=$root_dir/$install_prefix ;; # make absolute path
      esac
      ;;
    --build-dir) die "missing argument to $1 (see -h)" ;;
    --build-dir=*)
      build_dir=${1##*=}
      # Check if build_dir is an absolute path and if not, make it
      # absolute.
      case $build_dir in
        /*) ;;                               # absolute path
        *) build_dir=$root_dir/$build_dir ;; # make absolute path
      esac
      ;;
    --smt-switch-dir) die "missing argument to $1 (see -h)" ;;
    --smt-switch-dir=*)
      smt_switch_dir=${1##*=}
      # Check if this is an absolute path and if not, make it absolute.
      case $smt_switch_dir in
        /*) ;;                                         # absolute path
        *) smt_switch_dir=$root_dir/$smt_switch_dir ;; # make absolute path
      esac
      ;;
    --with-btor) with_boolector=ON ;;
    --with-msat) with_msat=ON ;;
    --with-yices2) with_yices2=ON ;;
    --with-z3) with_z3=ON ;;
    --with-msat-ic3ia) with_msat_ic3ia=ON ;;
    --with-coreir) with_coreir=ON ;;
    --with-coreir-extern) with_coreir_extern=ON ;;
    --debug)
      build_type=Debug
      ;;
    --docs) docs=ON ;;
    --python)
      python=ON
      ;;
    --static-lib)
      lib_type=STATIC
      ;;
    --static)
      static_exec=YES
      lib_type=STATIC
      ;;
    --with-profiling) with_profiling=ON ;;
    --no-system-gtest) system_gtest=OFF ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

[[ $lib_type == STATIC ]] && { [[ $with_coreir == ON ]] || [[ $with_coreir_extern == ON ]]; } &&
  die "CoreIR and static build are incompatible, must omit either '--static/--static-lib' or '--with-coreir/--with-coreir-extern'"

# Every option below is passed unconditionally, not just when it differs from
# its default, so that re-running configure.sh on an existing build
# directory always overwrites CMake's cache with the currently requested
# value (e.g. omitting --with-z3 on a later run correctly turns it back off,
# instead of leaving an earlier run's ON stuck in the cache). That in turn
# means we don't have to wipe the build directory to get a correct
# reconfigure.
cmake_opts=(
  -DCMAKE_BUILD_TYPE="$build_type"
  -DPONO_LIB_TYPE="$lib_type"
  -DPONO_STATIC_EXEC="$static_exec"
  -DSMT_SWITCH_DIR="$smt_switch_dir"
  -DCMAKE_INSTALL_PREFIX="$install_prefix"
  -DWITH_BOOLECTOR="$with_boolector"
  -DWITH_MSAT="$with_msat"
  -DWITH_YICES2="$with_yices2"
  -DWITH_Z3="$with_z3"
  -DWITH_MSAT_IC3IA="$with_msat_ic3ia"
  -DWITH_COREIR="$with_coreir"
  -DWITH_COREIR_EXTERN="$with_coreir_extern"
  -DBUILD_DOCS="$docs"
  -DBUILD_PYTHON_BINDINGS="$python"
  -DWITH_PROFILING="$with_profiling"
  -DSYSTEM_GTEST="$system_gtest"
)

echo "Running with cmake options: ${cmake_opts[*]}"

cmake -S "$root_dir" -B "$build_dir" "${cmake_opts[@]}" 2>&1
