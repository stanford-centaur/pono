#!/usr/bin/env bash
set -euo pipefail

smt_switch_version=e1369c7e078d2a1d72b3d6781285675ae321a10b

usage() {
  cat <<EOF
Usage: $0 [option] ...

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

-h, --help     display this message and exit
--with-btor    include Boolector (default: off)
--with-msat    include MathSAT which is under a custom non-BSD compliant license (default: off)
--with-yices2  include Yices2 which is under a custom non-BSD compliant license (default: off)
--with-z3      include Z3 (default: off)
--cvc5-home    use an already downloaded version of cvc5
--python       build python bindings (default: off)
EOF
  exit 0
}

die() {
  echo "*** $0: $*" >&2
  exit 1
}

# Navigate to source root
working_dir=$(pwd)
cd "$(dirname "${BASH_SOURCE[0]}")"
cd ..

# Parse options
conf_opts=()
with_boolector=false
with_z3=false
cvc5_home=default
incl_solver_str="bitwuzla and cvc5"
num_of_libs=3 # base, bitwuzla, cvc5
while (($# > 0)); do
  case $1 in
    -h | --help) usage ;;
    --with-msat)
      conf_opts+=(--msat --msat-home="$(pwd)/deps/mathsat")
      incl_solver_str="mathsat, $incl_solver_str"
      ((num_of_libs++))
      ;;
    --with-btor)
      with_boolector=true
      conf_opts+=(--btor)
      incl_solver_str="boolector, $incl_solver_str"
      ((num_of_libs++))
      ;;
    --with-yices2)
      conf_opts+=(--yices2)
      incl_solver_str="yices2, $incl_solver_str"
      ((num_of_libs++))
      ;;
    --with-z3)
      with_z3=true
      conf_opts+=(--z3)
      incl_solver_str="z3, $incl_solver_str"
      ((num_of_libs++))
      ;;
    --python)
      conf_opts+=(--python)
      ;;
    --cvc5-home) die "missing argument to $1 (see -h)" ;;
    --cvc5-home=*)
      cvc5_home="${1##*=}"
      # Check if cvc5_home is an absolute path and
      # if not, make it absolute.
      case $cvc5_home in
        /*) ;;
        *) cvc5_home="${working_dir}/$cvc5_home" ;;
      esac
      conf_opts+=(--cvc5-home="$cvc5_home")
      ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

# Download code if needed and check out correct version
mkdir -p deps
cd deps
if [[ ! -d smt-switch ]]; then
  git clone https://github.com/stanford-centaur/smt-switch
fi
cd smt-switch
git checkout -f "$smt_switch_version" || echo "warning: smt-switch folder is not a git repo"

# Build dependencies
./contrib/setup-bitwuzla.sh
if [[ $cvc5_home == default ]]; then
  ./contrib/setup-cvc5.sh
fi
if [[ $with_boolector == true ]]; then
  ./contrib/setup-boolector.sh
fi
if [[ $with_z3 == true ]]; then
  ./contrib/setup-z3.sh
fi

# Configure, build, test, and install smt-switch
if [[ -d build ]]; then
  echo "$(pwd)/build already exists, please remove it manually if you want to reconfigure smt-switch"
else
  ./configure.sh --prefix=local --static --smtlib-reader --bitwuzla --cvc5 "${conf_opts[@]}"
fi
cd build
cmake --build . -j
ctest -j
cmake --install .
cd ..

# Check that library files are there
lib_files=(local/lib/libsmt-switch*)
if ((${#lib_files[@]} == num_of_libs)); then
  echo "Smt-switch with $incl_solver_str was successfully installed to $(pwd)/local."
  echo "You may now build pono with: ./configure.sh && cd build && cmake --build ."
else
  echo "Building smt-switch failed."
  echo "You might be missing some dependencies."
  echo "Please see the github page for installation instructions: https://github.com/stanford-centaur/smt-switch"
  exit 1
fi
