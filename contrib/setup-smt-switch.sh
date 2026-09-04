#!/usr/bin/env bash
set -euo pipefail

smt_switch_version=460b5fdc8a048cdc2196de5bb27be2c2c9e6f90a # v1.1.4

# Keep the invocation directory to resolve relative path arguments against,
# since they would otherwise be taken relative to the smt-switch checkout that
# smt-switch's own configure.sh runs from.
working_dir=$(pwd)
root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd)"
cd "$root_dir"

usage() {
  cat <<EOF
Usage: $0 [<option> ...]

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

Each dependency location below suppresses the download and build of that
dependency. A --*-dir takes an install prefix, a --*-home a source tree.

-h, --help              display this message and exit

Dependency locations:
--bitwuzla-dir=STR      use an existing bitwuzla installation
--btor-home=STR         use an existing Boolector source tree
--cvc5-home=STR         use an existing cvc5 source tree
--msat-dir=STR          use an existing MathSAT installation (default: $root_dir/deps/mathsat)
--yices2-home=STR       use an existing Yices2 source tree
--z3-dir=STR            use an existing Z3 installation

Optional features / backends (default: off for all):
--python                build python bindings
--with-btor             include Boolector
--with-msat             include MathSAT which is under a custom non-BSD compliant license
--with-yices2           include Yices2 which is under a custom non-BSD compliant license
--with-z3               include Z3
EOF
  exit 0
}

die() {
  echo "*** $0: $*" >&2
  exit 1
}

abspath() {
  case $1 in
    /*) printf '%s\n' "$1" ;;
    *) printf '%s\n' "$working_dir/$1" ;;
  esac
}

# Parse options. This only records what was requested; the smt-switch options
# are assembled afterwards so that they do not depend on argument order.
bitwuzla_dir=
btor_home=
cvc5_home=
msat_dir=
yices2_home=
z3_dir=
with_boolector=false
with_msat=false
with_python=false
with_yices2=false
with_z3=false
while (($# > 0)); do
  case $1 in
    -h | --help) usage ;;

    # Dependency locations
    --bitwuzla-dir) die "missing argument to $1 (see -h)" ;;
    --bitwuzla-dir=*) bitwuzla_dir=$(abspath "${1#*=}") ;;
    --btor-home) die "missing argument to $1 (see -h)" ;;
    --btor-home=*) btor_home=$(abspath "${1#*=}") ;;
    --cvc5-home) die "missing argument to $1 (see -h)" ;;
    --cvc5-home=*) cvc5_home=$(abspath "${1#*=}") ;;
    --msat-dir) die "missing argument to $1 (see -h)" ;;
    --msat-dir=*) msat_dir=$(abspath "${1#*=}") ;;
    --yices2-home) die "missing argument to $1 (see -h)" ;;
    --yices2-home=*) yices2_home=$(abspath "${1#*=}") ;;
    --z3-dir) die "missing argument to $1 (see -h)" ;;
    --z3-dir=*) z3_dir=$(abspath "${1#*=}") ;;

    # Optional features / backends
    --python) with_python=true ;;
    --with-btor) with_boolector=true ;;
    --with-msat) with_msat=true ;;
    --with-yices2) with_yices2=true ;;
    --with-z3) with_z3=true ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

# A solver location on its own would make smt-switch's own configure.sh enable
# that backend implicitly, which would not match the libraries counted below.
if [[ -n $btor_home && $with_boolector == false ]]; then
  die "--btor-home requires --with-btor"
fi
if [[ -n $msat_dir && $with_msat == false ]]; then
  die "--msat-dir requires --with-msat"
fi
if [[ -n $yices2_home && $with_yices2 == false ]]; then
  die "--yices2-home requires --with-yices2"
fi
if [[ -n $z3_dir && $with_z3 == false ]]; then
  die "--z3-dir requires --with-z3"
fi

# Assemble the options for smt-switch's configure.sh, the solvers to report,
# and the options the user has to repeat when configuring pono itself.
conf_opts=()
pono_opts=()
solvers=(bitwuzla cvc5)
num_of_libs=3 # base, bitwuzla, cvc5
if [[ -n $bitwuzla_dir ]]; then
  conf_opts+=(--bitwuzla-dir="$bitwuzla_dir")
  pono_opts+=(--bitwuzla-dir="$bitwuzla_dir")
fi
if [[ -n $cvc5_home ]]; then
  conf_opts+=(--cvc5-home="$cvc5_home")
fi
if [[ $with_boolector == true ]]; then
  conf_opts+=(--btor)
  pono_opts+=(--with-btor)
  solvers+=(boolector)
  : $((num_of_libs++))
  if [[ -n $btor_home ]]; then
    conf_opts+=(--btor-home="$btor_home")
  fi
fi
if [[ $with_msat == true ]]; then
  conf_opts+=(--msat)
  pono_opts+=(--with-msat)
  solvers+=(mathsat)
  : $((num_of_libs++))
  if [[ -n $msat_dir ]]; then
    pono_opts+=(--msat-dir="$msat_dir")
  else
    # smt-switch looks in its own deps by default; use pono's, which is where
    # ci-scripts/setup-msat.sh unpacks the download.
    msat_dir=$root_dir/deps/mathsat
  fi
  conf_opts+=(--msat-home="$msat_dir")
fi
if [[ $with_yices2 == true ]]; then
  conf_opts+=(--yices2)
  pono_opts+=(--with-yices2)
  solvers+=(yices2)
  : $((num_of_libs++))
  if [[ -n $yices2_home ]]; then
    conf_opts+=(--yices2-home="$yices2_home")
  fi
fi
if [[ $with_z3 == true ]]; then
  conf_opts+=(--z3)
  pono_opts+=(--with-z3)
  solvers+=(z3)
  : $((num_of_libs++))
  if [[ -n $z3_dir ]]; then
    # smt-switch still spells --z3-dir --z3-install-dir.
    conf_opts+=(--z3-install-dir="$z3_dir")
    pono_opts+=(--z3-dir="$z3_dir")
  fi
fi
if [[ $with_python == true ]]; then
  conf_opts+=(--python)
  pono_opts+=(--python)
fi

# Download code if needed and check out correct version
mkdir -p deps
cd deps
if [[ ! -d smt-switch ]]; then
  git clone https://github.com/stanford-centaur/smt-switch
fi
cd smt-switch
git checkout -f "$smt_switch_version" ||
  echo "warning: smt-switch folder is not a git repo"

# Build dependencies, skipping the ones we were given a location for
if [[ -z $bitwuzla_dir ]]; then
  ./contrib/setup-bitwuzla.sh
fi
if [[ -z $cvc5_home ]]; then
  ./contrib/setup-cvc5.sh
fi
if [[ $with_boolector == true && -z $btor_home ]]; then
  ./contrib/setup-boolector.sh
fi
if [[ $with_z3 == true && -z $z3_dir ]]; then
  ./contrib/setup-z3.sh
fi

# Configure, build, test, and install smt-switch
if [[ -d build ]]; then
  echo "$(pwd)/build already exists, please remove it manually if you" \
    "want to reconfigure smt-switch"
else
  ./configure.sh --prefix=local --static --smtlib-reader --bitwuzla --cvc5 \
    --bison-dir=../bison-install ${conf_opts[@]+"${conf_opts[@]}"}
fi
cd build
cmake --build . -j
cmake --install .
cd ..

# Check that library files are there
lib_files=(local/lib/libsmt-switch*)
if ((${#lib_files[@]} == num_of_libs)); then
  solver_list=$(
    IFS=,
    echo "${solvers[*]}"
  )
  configure_cmd=./configure.sh
  if ((${#pono_opts[@]} > 0)); then
    configure_cmd="$configure_cmd ${pono_opts[*]}"
  fi
  echo "Smt-switch with ${solver_list//,/, } was successfully installed to" \
    "$(pwd)/local."
  echo "You may now build pono with: $configure_cmd && cd build && cmake" \
    "--build ."
else
  echo "Building smt-switch failed."
  echo "You might be missing some dependencies."
  echo "Please see the github page for installation instructions:" \
    "https://github.com/stanford-centaur/smt-switch"
  exit 1
fi
