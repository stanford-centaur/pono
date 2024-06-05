#!/usr/bin/env python3
from __future__ import annotations

import argparse
import enum
import os
import signal
import subprocess
import sys
import time


class ReturnCode(enum.Enum):
    """Return codes corresponding to enum values in pono/core/proverresult.h."""

    SEGFAULT = -signal.SIGSEGV
    SAT = 0
    UNSAT = 1
    ERROR = 2
    UNKNOWN = 255


SOLVED_RETURN_CODES = {ReturnCode.SAT.value, ReturnCode.UNSAT.value}


def main():
    parser = argparse.ArgumentParser(description="Run multiple engines in parallel")
    parser.add_argument("btor_file")
    parser.add_argument(
        "-k", "--bound", default=1000, type=int, help="The maximum bound to unroll to"
    )
    parser.add_argument("-s", "--smt-solver", default="bzla", help="SMT solver to use")
    parser.add_argument(
        "-c", "--csv-summary", action="store_true", help="print csv summary to stderr"
    )
    args = parser.parse_args()

    def get_command(arguments: list[str]) -> list[str]:
        command = ["pono", "-k", str(args.bound), *arguments]
        if "interp" in arguments or "ic3ia" in arguments:
            # Interpolation requires MathSat
            command.extend(["--smt-solver", "msat"])
        else:
            command.extend(["--smt-solver", args.smt_solver])
            if args.smt_solver == "btor" and "--ceg-bv-arith" in command:
                # BV UF abstraction doesn't work with plain Boolector
                command.append("--logging-smt-solver")
        command.append(args.btor_file)
        return command

    pono_arguments = {
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
        "K-Induction": ["-e", "ind"],
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

    processes: dict[str, tuple[subprocess.Popen[str], float]] = {}

    def terminate_all():
        for process, _ in processes.values():
            if process.poll() is None:
                process.terminate()

    def handle_signal(signum, frame):
        # send signal received to subprocesses
        terminate_all()
        sys.exit(ReturnCode.UNKNOWN.value)

    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)

    for name, pono_args in pono_arguments.items():
        cmd = get_command(pono_args)
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=open(os.devnull, "w"), text=True
        )
        processes[name] = proc, time.time()

    retcode = ReturnCode.UNKNOWN.value

    while processes:
        for name, (process, start_time) in list(processes.items()):
            end_time = time.time()
            if process.poll() is not None:
                del processes[name]
                if process.returncode == ReturnCode.UNKNOWN.value:
                    continue
                if args.csv_summary:
                    try:
                        result = ReturnCode(process.returncode).name.lower()
                    except ValueError:
                        result = f"error({process.returncode})"
                    duration = end_time - start_time
                    cmd = " ".join(get_command(pono_arguments[name]))
                    print(
                        result, f"{duration:.1f}", name, cmd, sep=",", file=sys.stderr
                    )
                if process.returncode in SOLVED_RETURN_CODES:
                    retcode = process.returncode
                    assert process.stdout
                    print(process.stdout.read())
                    break
        else:
            time.sleep(0.001)
            continue
        terminate_all()
        break

    return retcode


if __name__ == "__main__":
    sys.exit(main())
