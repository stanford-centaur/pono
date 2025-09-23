#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import csv
import dataclasses
import enum
import logging
import pathlib
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from typing import cast


class ReturnCode(enum.Enum):
    """Return codes corresponding to enum values in pono/core/proverresult.h."""

    SAT = 0
    UNSAT = 1
    ERROR = 2
    UNKNOWN = 255  # originally -1, but negative POSIX return codes wrap around


SOLVED_RETURN_CODES = {ReturnCode.SAT.value, ReturnCode.UNSAT.value}


@dataclasses.dataclass
class Summary:
    """CSV summary record for a single pono run."""

    result: str
    runtime: str
    engine: str
    command: str


SUMMARY_FIELDS = [field.name for field in dataclasses.fields(Summary)]


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


def summarize(file: str, engine: str, returncode: int, runtime: float, cmd: list[str]):
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
    summary = Summary(
        result=result,
        runtime=f"{runtime:.1f}",
        engine=engine,
        command=" ".join(cmd),
    )
    with open(file, "a") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=SUMMARY_FIELDS)
        writer.writerow(dataclasses.asdict(summary))


def clean_up(
    processes: dict[str, subprocess.Popen[str]],
    witnesses: dict[str, pathlib.Path],
    verbose: bool,
):
    for name, process in processes.items():
        if process.poll() is None:  # process has not finished yet
            process.terminate()
        if name in witnesses:
            witnesses[name].unlink(missing_ok=True)
        if verbose and process.stderr and (stderr := process.stderr.read()):
            logger = logging.getLogger(name)
            for line in stderr.splitlines():
                logger.warning(line)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run multiple engines in parallel")
    parser.add_argument("btor_file", help="input benchmark in BTOR2 format")
    parser.add_argument("witness_file", nargs="?", help="file to store the witness")
    parser.add_argument("-k", "--bound", default=1000, type=int, help="unrolling bound")
    parser.add_argument("-v", "--verbose", action="store_true", help="echo stderr")
    parser.add_argument("-s", "--summarize", metavar="FILE", help="save csv summary")
    args = parser.parse_args()

    # Create summary file when needed, truncating if it exists.
    if args.summarize:
        with open(args.summarize, "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=SUMMARY_FIELDS)
            writer.writeheader()

    # Configure logging.
    logging.basicConfig(format="{name}: {message}", style="{")
    logger = logging.getLogger(parser.prog)

    # Try current directory if pono is not on PATH.
    executable = shutil.which("pono") or pathlib.Path(__file__).parent / "pono"

    processes: dict[str, subprocess.Popen[str]] = {}
    start_times: dict[str, float] = {}
    witnesses: dict[str, pathlib.Path] = {}
    atexit.register(clean_up, processes, witnesses, args.verbose)

    # Launch each portfolio solver as a subprocess.
    for name, options in ENGINE_OPTIONS.items():
        cmd = [executable, "-k", str(args.bound), *options]
        if args.witness_file:
            witness_file = tempfile.NamedTemporaryFile(delete=False)
            witness_file.close()
            witnesses[name] = pathlib.Path(witness_file.name)
            cmd.extend(["--dump-btor2-witness", witness_file.name])
        cmd.append(args.btor_file)
        stderr = subprocess.PIPE if args.verbose else subprocess.DEVNULL
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr, text=True)
        start_times[name] = time.time()
        processes[name] = proc

    # Wait for at least one process to terminate successfully.
    while processes:
        for name, process in processes.items():
            end_time = time.time()
            if process.poll() is not None:  # process has finished
                if args.summarize:
                    runtime = end_time - start_times[name]
                    cmd = cast(list[str], process.args)
                    summarize(args.summarize, name, process.returncode, runtime, cmd)
                if process.returncode in SOLVED_RETURN_CODES:
                    if process.stdout is None:
                        logger.warning(f"{name} has no stdout")
                        print(ReturnCode(process.returncode).name.lower())
                    else:
                        print(process.stdout.read())
                    if args.witness_file:
                        witnesses[name].rename(args.witness_file)
                    return process.returncode
                else:
                    del processes[name]
                    clean_up({name: process}, witnesses, args.verbose)
                    break

    return ReturnCode.UNKNOWN.value


if __name__ == "__main__":
    sys.exit(main())
