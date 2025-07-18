#!/bin/sh

# Syntax and structure borrowed from cvc5's configure.sh script

CONF_FILE=Makefile.conf

usage () {
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

die () {
    echo "*** $0: $*" 1>&2
    exit 1
}

build_dir=build
smt_switch_dir=default
install_prefix=default
build_type=default
with_boolector=default
with_msat=default
with_yices2=default
with_msat_ic3ia=default
with_coreir=default
with_coreir_extern=default
debug=default
docs=default
python=default
lib_type=SHARED
static_exec=NO
with_profiling=default
system_gtest=default

buildtype=Release

while [ $# -gt 0 ]
do
    case $1 in
        -h|--help) usage;;
        --prefix) die "missing argument to $1 (see -h)" ;;
        --prefix=*)
            install_prefix=${1##*=}
            # Check if install_prefix is an absolute path and if not, make it
            # absolute.
            case $install_prefix in
                /*) ;;                                      # absolute path
                *) install_prefix=$(pwd)/$install_prefix ;; # make absolute path
            esac
            ;;
        --build-dir) die "missing argument to $1 (see -h)" ;;
        --build-dir=*)
            build_dir=${1##*=}
            # Check if build_dir is an absolute path and if not, make it
            # absolute.
            case $build_dir in
                /*) ;;                                      # absolute path
                *) build_dir=$(pwd)/$build_dir ;; # make absolute path
            esac
            ;;
        --smt-switch-dir) die "missing argument to $1 (see -h)" ;;
        --smt-switch-dir=*)
            smt_switch_dir=${1##*=}
            # Check if this is an absolute path and if not, make it absolute.
            case $smt_switch_dir in
                /*) ;;                                      # absolute path
                *) smt_switch_dir=$(pwd)/$smt_switch_dir ;; # make absolute path
            esac
            ;;
        --with-btor) with_boolector=ON;;
        --with-msat) with_msat=ON;;
        --with-yices2) with_yices2=ON;;
        --with-msat-ic3ia) with_msat_ic3ia=ON;;
        --with-coreir) with_coreir=ON;;
        --with-coreir-extern) with_coreir_extern=ON;;
        --debug)
            debug=yes;
            buildtype=Debug
            ;;
        --docs) docs=yes;;
        --python)
            python=yes
            ;;
        --static-lib)
            lib_type=STATIC
            ;;
        --static)
            static_exec=YES;
            lib_type=STATIC;
            ;;
        --with-profiling) with_profiling=ON;;
        --no-system-gtest) system_gtest=no;;
        *) die "unexpected argument: $1";;
    esac
    shift
done

[ $lib_type = STATIC ] && [ $with_coreir = ON -o $with_coreir_extern = ON ] && \
    die "CoreIR and static build are incompatible, must omit either '--static/--static-lib' or '--with-coreir/--with-coreir-extern'"

cmake_opts="-DCMAKE_BUILD_TYPE=$buildtype -DPONO_LIB_TYPE=${lib_type} -DPONO_STATIC_EXEC=${static_exec}"

[ $smt_switch_dir != default ] \
    && cmake_opts="$cmake_opts -DSMT_SWITCH_DIR=${smt_switch_dir}"

[ $install_prefix != default ] \
    && cmake_opts="$cmake_opts -DCMAKE_INSTALL_PREFIX=$install_prefix"

[ $with_boolector != default ] \
    && cmake_opts="$cmake_opts -DWITH_BOOLECTOR=$with_boolector"

[ $with_msat != default ] \
    && cmake_opts="$cmake_opts -DWITH_MSAT=$with_msat"

[ $with_yices2 != default ] \
    && cmake_opts="$cmake_opts -DWITH_YICES2=$with_yices2"

[ $with_msat_ic3ia != default ] \
    && cmake_opts="$cmake_opts -DWITH_MSAT_IC3IA=$with_msat_ic3ia"

[ $with_coreir != default ] \
    && cmake_opts="$cmake_opts -DWITH_COREIR=$with_coreir"

[ $with_coreir_extern != default ] \
    && cmake_opts="$cmake_opts -DWITH_COREIR_EXTERN=$with_coreir_extern"

[ $docs != default ] \
    && cmake_opts="$cmake_opts -DBUILD_DOCS=ON"

[ $python != default ] \
    && cmake_opts="$cmake_opts -DBUILD_PYTHON_BINDINGS=ON"

[ $with_profiling != default ] \
    && cmake_opts="$cmake_opts -DWITH_PROFILING=$with_profiling"

[ $system_gtest != default ] \
    && cmake_opts="$cmake_opts -DSYSTEM_GTEST=$system_gtest"

root_dir=$(pwd)

[ -e "$build_dir" ] && rm -r "$build_dir"

mkdir -p "$build_dir"
cd "$build_dir" || exit 1

[ -e CMakeCache.txt ] && rm CMakeCache.txt

echo "Running with cmake options: $cmake_opts"

cmake "$root_dir" $cmake_opts 2>&1
