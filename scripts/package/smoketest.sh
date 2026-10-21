#!/usr/bin/env bash
# Smoke test of the pono binary beside this script, which SV-COMP requires of
# a tool archive: prints the version and the usage of pono, then runs it on
# the MoXI tasks in examples/, each named after its expected verdict. Exits
# with status 0 if and only if pono reports every expected verdict.
#
# Set PONO to test another binary, e.g. PONO=build/pono in the source tree.
set -euo pipefail

dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
pono=${PONO:-$dir/pono}
# IC3 with implicit predicate abstraction, with MathSAT as the solver and the
# interpolator. Bounded model checking, the default engine, cannot prove a
# query unsatisfiable, and Bitwuzla, the default solver, does not support
# integers.
options=(-e ic3ia --smt-solver msat --smt-interpolator msat)

# pono exits with status 2 after printing its version or usage.
"$pono" --version || [[ $? -eq 2 ]]
"$pono" --help || [[ $? -eq 2 ]]

failed=0
for expected in sat unsat; do
  echo
  echo "\$ pono ${options[*]} examples/$expected.moxi"
  # pono exits with status 0 for sat and 1 for unsat, so the verdict is the
  # last line of its output that starts with one, as BenchExec reads it.
  output=$("$pono" "${options[@]}" "$dir/examples/$expected.moxi") || true
  echo "$output"
  verdict=$(grep -E '^(sat|unsat|unknown)' <<<"$output" | tail -n 1) || true
  if [[ $verdict == "$expected" ]]; then
    echo "OK: $expected, as expected"
  else
    echo "FAILED: expected $expected, got ${verdict:-no verdict}" >&2
    failed=1
  fi
done
exit "$failed"
