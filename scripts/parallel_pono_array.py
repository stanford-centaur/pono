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
    "bmc": [
        "-e",
        "bmc",
        "--static-coi",
        "--bmc-allow-non-minimal-cex",
        "--bmc-bound-start",
        "5",
        "--bmc-exponential-step",
    ],
    "ind": [
        "-e",
        "ind",
        "--static-coi",
    ],
    "ind.ceg_bv_arith": [
        "-e",
        "ind",
        "--static-coi",
        "--ceg-bv-arith",
    ],
    "ind.ceg_bv_arith.bound_step.11": [
        "-e",
        "ind",
        "--static-coi",
        "--ceg-bv-arith",
        "--kind-no-simple-path-check",
        "--kind-bound-step",
        "11",
        "--kind-one-time-base-check",
    ],
    "interp": [
        "-e",
        "interp",
        "--smt-interpolator",
        "bzla",
        "--static-coi",
        "--promote-inputvars",
    ],
    "interp.backward": [
        "-e",
        "interp",
        "--smt-interpolator",
        "bzla",
        "--static-coi",
        "--promote-inputvars",
        "--interp-backward",
    ],
    "interp.eager_unroll": [
        "-e",
        "interp",
        "--smt-interpolator",
        "bzla",
        "--static-coi",
        "--promote-inputvars",
        "--interp-eager-unroll",
    ],
    "interp.ceg_bv_arith": [
        "-e",
        "interp",
        "--smt-interpolator",
        "bzla",
        "--static-coi",
        "--promote-inputvars",
        "--ceg-bv-arith",
    ],
    "interp.ceg_bv_arith.backward": [
        "-e",
        "interp",
        "--smt-interpolator",
        "bzla",
        "--static-coi",
        "--promote-inputvars",
        "--ceg-bv-arith",
        "--interp-backward",
    ],
    "ismc": [
        "-e",
        "ismc",
        "--static-coi",
    ],
    "ismc.ceg_bv_arith": [
        "-e",
        "ismc",
        "--static-coi",
        "--ceg-bv-arith",
        "--promote-inputvars",
    ],
    "ic3ia": [
        "-e",
        "ic3ia",
        "--smt-solver",
        "msat",
        "--pseudo-init-prop",
        "--ceg-prophecy-arrays",
    ],
    "ic3ia.ceg_bv_arith": [
        "-e",
        "ic3ia",
        "--smt-solver",
        "msat",
        "--pseudo-init-prop",
        "--ceg-prophecy-arrays",
        "--ceg-bv-arith",
    ],
}


def summarize(
    file: pathlib.Path,
    engine: str,
    returncode: int,
    runtime: float,
    cmd: list[str],
) -> None:
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
    with file.open("a") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=SUMMARY_FIELDS)
        writer.writerow(dataclasses.asdict(summary))


def clean_up(
    processes: dict[str, subprocess.Popen[str]],
    witnesses: dict[str, pathlib.Path],
    *,
    verbose: bool,
) -> None:
    for name, process in processes.items():
        if process.poll() is None:  # process has not finished yet
            process.terminate()
        if name in witnesses:
            witnesses[name].unlink(missing_ok=True)
        if verbose and process.stderr and (stderr := process.stderr.read()):
            logger = logging.getLogger(name)
            for line in stderr.splitlines():
                logger.warning(line)


def find_file(name: str) -> pathlib.Path:
    """Return the path to a file in PATH or the script's directory."""
    fullname = shutil.which(name)
    if fullname is None:
        path = pathlib.Path(__file__).parent / name
    else:
        path = pathlib.Path(fullname)
    if not path.is_file():
        msg = f"file {name} not  in script directory or PATH"
        raise FileNotFoundError(msg)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run multiple engines in parallel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("btor_file", help="input benchmark in BTOR2 format")
    parser.add_argument("witness_file", nargs="?", help="file to store the witness in")
    parser.add_argument("-b", "--binary", default="pono", help="name of pono binary")
    parser.add_argument("-k", "--bound", default=2**20, type=int, help="check until")
    parser.add_argument("-v", "--verbose", action="store_true", help="echo stderr")
    parser.add_argument(
        "-s",
        "--summarize",
        metavar="FILE",
        type=pathlib.Path,
        help="save a csv summary to the specified file",
    )
    args = parser.parse_args()

    # Create summary file when needed, truncating if it exists.
    if args.summarize:
        with args.summarize.open("w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=SUMMARY_FIELDS)
            writer.writeheader()

    # Configure logging.
    logging.basicConfig(format="{name}: {message}", style="{")
    logger = logging.getLogger(parser.prog)

    processes: dict[str, subprocess.Popen[str]] = {}
    start_times: dict[str, float] = {}
    witnesses: dict[str, pathlib.Path] = {}
    atexit.register(clean_up, processes, witnesses, verbose=args.verbose)

    # Launch each portfolio solver as a subprocess.
    executable = find_file(args.binary)
    for name, options in ENGINE_OPTIONS.items():
        cmd = [str(executable), "-k", str(args.bound), *options]
        if args.witness_file:
            with tempfile.NamedTemporaryFile(delete=False) as witness_file:
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
                    cmd = cast("list[str]", process.args)
                    summarize(args.summarize, name, process.returncode, runtime, cmd)
                if process.returncode in SOLVED_RETURN_CODES:
                    if process.stdout is None:
                        logger.warning("%s has no stdout", name)
                        print(ReturnCode(process.returncode).name.lower())
                    else:
                        print(process.stdout.read())
                    if args.witness_file:
                        shutil.move(witnesses[name], args.witness_file)
                    return process.returncode
                del processes[name]
                clean_up({name: process}, witnesses, verbose=args.verbose)
                break

    return ReturnCode.UNKNOWN.value


if __name__ == "__main__":
    sys.exit(main())
