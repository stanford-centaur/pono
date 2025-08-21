#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import enum
import shutil
import signal
import subprocess
import sys
import time

ProcessMap = dict[str, tuple[list[str], subprocess.Popen[str], float]]


class ReturnCode(enum.Enum):
    """Return codes corresponding to enum values in pono/core/proverresult.h."""

    SAT = 0
    UNSAT = 1
    ERROR = 2
    UNKNOWN = 255  # originally -1, but negative POSIX return codes wrap around


SOLVED_RETURN_CODES = {ReturnCode.SAT.value, ReturnCode.UNSAT.value}

ENGINE_OPTIONS = {
    "BMC": [
        "-e",
        "bmc",
        "--static-coi",
        "--bmc-bound-start",
        "0",
        "--bmc-bound-step",
        "11",
        "--bmc-single-bad-state",
    ],
    "K-Induction": ["-e", "ind", "--static-coi"],
    "Interpolant-based": ["-e", "interp"],
    "MBIC3": ["-e", "mbic3", "--static-coi"],
    "IC3Bits": ["-e", "ic3bits"],
    "IC3IA": ["-e", "ic3ia", "--pseudo-init-prop"],
    "IC3IA-UF": ["-e", "ic3ia", "--pseudo-init-prop", "--ceg-bv-arith"],
    "IC3IA-NoUCG": ["-e", "ic3ia", "--pseudo-init-prop", "--no-ic3-unsatcore-gen"],
    "IC3IA-FTS": ["-e", "ic3ia", "--static-coi"],
    "IC3SA": ["-e", "ic3sa", "--static-coi"],
    "IC3SA-UF": ["-e", "ic3sa", "--static-coi", "--ceg-bv-arith"],
    "SyGuS-PDR": ["-e", "sygus-pdr", "--promote-inputvars"],
    "SyGuS-PDR-2": [
        "-e",
        "sygus-pdr",
        "--promote-inputvars",
        "--sygus-term-mode",
        "2",
    ],
}


def print_summary(engine: str, returncode: int, runtime: float, cmd: list[str]):
    if returncode < 0:
        signum = -returncode
        try:
            result = signal.Signals(signum).name.lower()
        except ValueError:
            result = f"signal {signum}"
    else:
        try:
            result = ReturnCode(returncode).name.lower()
        except ValueError:
            result = f"error {returncode}"
    print(result, f"{runtime:.1f}", engine, " ".join(cmd), sep=",", file=sys.stderr)


def clean_up(processes: ProcessMap, verbose: bool):
    for name, (_, process, _) in processes.items():
        if process.poll() is None:
            process.terminate()
        if verbose and process.stderr and (stderr := process.stderr.read()):
            for line in stderr.splitlines():
                print(f"{name}:", line, file=sys.stderr)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run multiple engines in parallel")
    parser.add_argument("btor_file", help="input benchmark in BTOR2 format")
    parser.add_argument("witness_file", nargs="?", help="file to store the witness")
    parser.add_argument("-k", "--bound", default=1000, type=int, help="unrolling bound")
    parser.add_argument("-s", "--smt-solver", help="main SMT solver")
    parser.add_argument("-v", "--verbose", action="store_true", help="echo stderr")
    parser.add_argument("--summarize", action="store_true", help="print summary report")
    args = parser.parse_args()

    processes: dict[str, tuple[list[str], subprocess.Popen[str], float]] = {}
    atexit.register(clean_up, processes, args.verbose)

    # Try current directory if pono is not on PATH.
    executable = shutil.which("pono") or "./pono"

    for name, options in ENGINE_OPTIONS.items():
        cmd = [executable, "-k", str(args.bound), *options]
        if args.witness_file:
            cmd.extend(["--dump-btor2-witness", args.witness_file])
        if args.smt_solver:
            cmd.extend(["--smt-solver", args.smt_solver])
        if args.smt_solver == "btor" and "--ceg-bv-arith" in cmd:
            # BV UF abstraction doesn't work with plain Boolector
            cmd.append("--logging-smt-solver")
        cmd.append(args.btor_file)
        stderr = subprocess.PIPE if args.verbose else subprocess.DEVNULL
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr, text=True)
        processes[name] = cmd, proc, time.time()

    while processes:
        for name, (cmd, process, start_time) in list(processes.items()):
            end_time = time.time()
            if process.poll() is not None:
                del processes[name]
                clean_up({name: ([], process, 0)}, args.verbose)
                if args.summarize:
                    print_summary(name, process.returncode, end_time - start_time, cmd)
                if process.returncode in SOLVED_RETURN_CODES:
                    if process.stdout is None:
                        print(f"{parser.prog}:", name, "has no stdout", file=sys.stderr)
                        print(ReturnCode(process.returncode).name.lower())
                    else:
                        print(process.stdout.read())
                    return process.returncode

    return ReturnCode.UNKNOWN.value


if __name__ == "__main__":
    sys.exit(main())
