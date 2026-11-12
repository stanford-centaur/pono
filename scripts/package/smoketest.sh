#!/usr/bin/env bash
# Smoke test of the package beside this script, which SV-COMP requires of a
# tool archive: prints the version and the usage of pono, the portfolio
# script, then runs its portfolio for SV-COMP on the MoXI tasks in examples/,
# each named after its expected verdict, as SV-COMP runs it. Exits with
# status 0 if and only if pono reports every expected verdict.
#
# Set PONO_BIN to run the portfolio with another binary than pono-bin, e.g.
# PONO_BIN=build/pono scripts/package/smoketest.sh in the source tree.
set -euo pipefail

dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
pono=$dir/pono
binary=()
if [[ -n ${PONO_BIN:-} ]]; then
  binary=(--pono-bin "$PONO_BIN")
fi

# pono exits with status 2 after printing the version of pono-bin, as pono-bin
# does.
"$pono" "${binary[@]}" --version || [[ $? -eq 2 ]]
"$pono" --help

failed=0
for expected in sat unsat; do
  echo
  echo "\$ pono --svcomp examples/$expected.moxi"
  # pono exits with status 0 for sat and 1 for unsat, so the verdict is the
  # last line of its output that starts with one, as BenchExec reads it.
  output=$("$pono" --svcomp "${binary[@]}" "$dir/examples/$expected.moxi") ||
    true
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
