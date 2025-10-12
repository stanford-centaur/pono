#!/usr/bin/env bash
set -euo pipefail

script_dir="$(dirname "${BASH_SOURCE[0]}")"
deps_dir="$(dirname "${script_dir}")/deps"
smt_switch_dir="{deps_dir}/smt-switch"

smt_switch_version=v1.0.1

usage() {
  cat <<EOF
Usage: $0 [option] ...

Sets up the smt-switch API for interfacing with SMT solvers through a C++ API.

-h, --help     display this message and exit
--with-btor    include Boolector (default: off)
--with-msat    include MathSAT which is under a custom non-BSD compliant license (default: off)
--with-yices2  include Yices2 which is under a custom non-BSD compliant license (default: off)
--cvc5-home    use an already downloaded version of cvc5
--python       build python bindings (default: off)
EOF
  exit 0
}

die() {
  echo "*** $0: $*" >&2
  exit 1
}

with_boolector=default
conf_opts=()
cvc5_home=default
incl_solver_str="bitwuzla and cvc5"
num_of_libs=3 # base, bitwuzla, cvc5
while (($# > 0)); do
  case $1 in
    -h | --help) usage ;;
    --with-msat)
      conf_opts+=(--msat --msat-home=../mathsat)
      incl_solver_str="mathsat, ${incl_solver_str}"
      ((num_of_libs++))
      ;;
    --with-btor)
      with_boolector=ON
      conf_opts+=(--btor)
      incl_solver_str="boolector, ${incl_solver_str}"
      ((num_of_libs++))
      ;;
    --with-yices2)
      conf_opts+=(--yices2)
      incl_solver_str="yices2, ${incl_solver_str}"
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
        *) cvc5_home="$(pwd)/${cvc5_home}" ;;
      esac
      conf_opts+=(--cvc5-home="${cvc5_home}")
      ;;
    *) die "unexpected argument: $1" ;;
  esac
  shift
done

if [[ ! -d ${smt_switch_dir} ]]; then
  mkdir -p "${deps_dir}"
  cd "${deps_dir}"
  git clone https://github.com/stanford-centaur/smt-switch
  cd smt-switch
  git checkout -f "${smt_switch_version}"
  ./contrib/setup-bitwuzla.sh
  if [[ ${cvc5_home} == default ]]; then
    ./contrib/setup-cvc5.sh
  fi
  if [[ ${with_boolector} == ON ]]; then
    ./contrib/setup-boolector.sh
  fi
  ./configure.sh --bitwuzla --cvc5 "${conf_opts[@]}" --prefix=local --static --smtlib-reader
  cd build
  make -j
  make test
  make install
else
  echo "${smt_switch_dir} already exists. If you want to rebuild, please remove it manually."
fi

lib_files=("${smt_switch_dir}"/smt-switch/local/lib/libsmt-switch*)
if ((${#lib_files[@]} != num_of_libs)); then
  echo "It appears smt-switch with $incl_solver_str was successfully installed to ${smt_switch_dir}/local."
  echo "You may now build pono with: ./configure.sh && cd build && make"
else
  echo "Building smt-switch failed."
  echo "You might be missing some dependencies."
  echo "Please see the github page for installation instructions: https://github.com/stanford-centaur/smt-switch"
  exit 1
fi
