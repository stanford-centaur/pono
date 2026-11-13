#!/usr/bin/env bash
# Smoke test of the package beside this script, which SV-COMP requires of a
# tool archive: prints the version and the usage of pono, the portfolio
# script, then runs on each MoXI task in examples/, which is named after its
# expected verdict, its portfolio for SV-COMP, as SV-COMP runs it, and
# pono-bin with each SMT solver that the portfolio uses. Exits with status 0
# if and only if every run reports the expected verdict.
#
# Set PONO_BIN to test another binary than pono-bin, e.g.
# PONO_BIN=build/pono scripts/package/smoketest.sh in the source tree.
set -euo pipefail

dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
pono=$dir/pono
pono_bin=${PONO_BIN:-$dir/pono-bin}
binary=()
if [[ -n ${PONO_BIN:-} ]]; then
  binary=(--pono-bin "$PONO_BIN")
fi

# pono exits with status 2 after printing the version of pono-bin, as pono-bin
# does.
"$pono" "${binary[@]}" --version || [[ $? -eq 2 ]]
"$pono" --help

failed=0

# Runs the command $2... on the task with the expected verdict $1.
check() {
  local expected=$1 output verdict
  shift
  # pono and pono-bin exit with status 0 for sat and 1 for unsat, so the
  # verdict is the last line of the output that starts with one, as BenchExec
  # reads it.
  output=$("$@" "$dir/examples/$expected.moxi") || true
  echo "$output"
  verdict=$(grep -E '^(sat|unsat|unknown)' <<<"$output" | tail -n 1) || true
  if [[ $verdict == "$expected" ]]; then
    echo "OK: $expected, as expected"
  else
    echo "FAILED: expected $expected, got ${verdict:-no verdict}" >&2
    failed=1
  fi
}

for expected in sat unsat; do
  echo
  echo "\$ pono --svcomp examples/$expected.moxi"
  check "$expected" "$pono" --svcomp "${binary[@]}"
  for solver in bzla cvc5 msat z3; do
    options=(--engine ind --smt-solver "$solver")
    echo
    echo "\$ pono-bin ${options[*]} examples/$expected.moxi"
    check "$expected" "$pono_bin" "${options[@]}"
  done
done
exit "$failed"
