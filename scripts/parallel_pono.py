#!/usr/bin/env python3
from __future__ import annotations

import argparse
import atexit
import csv
import dataclasses
import enum
import logging
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from logging import Logger
from pathlib import Path
from typing import NamedTuple, cast


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


def summarize(
    file: Path,
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


def find_executable(name: str) -> Path:
    """Return the path to a file in PATH or the script's directory."""
    fullname = shutil.which(name)
    path = Path(__file__).parent / name if fullname is None else Path(fullname)
    if not path.is_file():
        msg = f"file '{name}' is not in script directory or PATH"
        raise FileNotFoundError(msg)
    if not os.access(path, os.X_OK):
        msg = f"{path} is not executable"
        raise PermissionError(msg)
    return path


class ProverDesc(NamedTuple):
    name: str
    process: subprocess.Popen
    start_time: float
    witness_file: Path | None = None


def clean_up_prover(desc: ProverDesc, *, verbose: bool) -> None:
    if desc.process.poll() is None:  # process has not finished yet
        desc.process.terminate()
    if desc.witness_file:
        desc.witness_file.unlink(missing_ok=True)
    if verbose and desc.process.stderr:
        stderr = desc.process.stderr.read()
        logger = logging.getLogger(desc.name)
        for line in stderr.splitlines():
            logger.warning(line)


class ProverManager:
    def __init__(self, binary_name: str, logger: Logger, *, verbose: bool) -> None:
        self.provers: set[ProverDesc] = set()
        self.executable = str(find_executable(binary_name))
        self.logger = logger
        self.verbose = verbose
        atexit.register(self._clean_up)

    def _clean_up(self) -> None:
        for prover in self.provers:
            clean_up_prover(prover, verbose=self.verbose)

    def start_provers(self, btor_file: str, bound: int) -> None:
        """Launch each portfolio solver as a subprocess."""
        for name, options in ENGINE_OPTIONS.items():
            cmd = [self.executable, "-k", str(bound), *options, btor_file]
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE if self.verbose else subprocess.DEVNULL,
                text=True,
            )
            self.provers.add(ProverDesc(name, proc, time.time()))

    def start_provers_with_witness(self, btor_file: str, bound: int) -> None:
        """Launch each portfolio solver as a subprocess."""
        for name, options in ENGINE_OPTIONS.items():
            witness_file = tempfile.NamedTemporaryFile(delete=False)  # noqa: SIM115
            cmd = [
                self.executable,
                "-k",
                str(bound),
                *options,
                "--dump-btor2-witness",
                witness_file.name,
                btor_file,
            ]
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE if self.verbose else subprocess.DEVNULL,
                text=True,
            )
            self.provers.add(
                ProverDesc(name, proc, time.time(), Path(witness_file.name))
            )

    def wait_for_provers(self, summary_file: Path | None) -> tuple[int, Path | None]:
        while self.provers:
            for prover in self.provers:
                end_time = time.time()
                if prover.process.poll() is not None:  # process has finished
                    if summary_file:
                        runtime = end_time - prover.start_time
                        cmd = cast("list[str]", prover.process.args)
                        summarize(
                            summary_file,
                            prover.name,
                            prover.process.returncode,
                            runtime,
                            cmd,
                        )
                    if prover.process.returncode not in SOLVED_RETURN_CODES:
                        self.provers.discard(prover)
                        clean_up_prover(prover, verbose=self.verbose)
                        break
                    if prover.process.stdout is None:
                        self.logger.warning("%s has no stdout", prover.name)
                        print(ReturnCode(prover.process.returncode).name.lower())
                    else:
                        print(prover.process.stdout.read())
                    return prover.process.returncode, prover.witness_file
        return ReturnCode.UNKNOWN.value, None


def get_args() -> tuple[argparse.Namespace, str]:
    parser = argparse.ArgumentParser(
        description="Run multiple engines in parallel",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("btor_file", help="input benchmark in BTOR2 format")
    parser.add_argument("witness_file", nargs="?", help="file to store the witness in")
    parser.add_argument(
        "-b",
        "--binary",
        default="pono",
        help="name of pono binary (file must be in PATH or script dir)",
    )
    parser.add_argument(
        "-k", "--bound", default=2**20, type=int, help="bound to check until"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="echo stderr")
    parser.add_argument(
        "-s",
        "--summarize",
        metavar="FILE",
        type=Path,
        help="save a csv summary to the specified file",
    )
    return parser.parse_args(), parser.prog


def main() -> int:
    args, progname = get_args()
    # Create summary file when needed, truncating if it exists.
    if args.summarize:
        with args.summarize.open("w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=SUMMARY_FIELDS)
            writer.writeheader()

    # Configure logging.
    logging.basicConfig(format="{name}: {message}", style="{")
    logger = logging.getLogger(progname)

    prover_manager = ProverManager(args.binary, logger, verbose=args.verbose)
    if args.witness_file:
        prover_manager.start_provers_with_witness(args.btor_file, args.bound)
    else:
        prover_manager.start_provers(args.btor_file, args.bound)
    return_code, witness_file = prover_manager.wait_for_provers(args.summarize)

    if args.witness_file:
        if not witness_file:
            logger.warning("witness requested but no witness file present")
        else:
            shutil.move(witness_file, args.witness_file)

    return return_code


if __name__ == "__main__":
    sys.exit(main())
