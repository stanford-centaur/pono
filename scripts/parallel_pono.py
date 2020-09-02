#!/usr/bin/env python3

import argparse
import signal
import subprocess
import sys
import time
import os

## Non-blocking reads for subprocess
## https://stackoverflow.com/questions/375427/non-blocking-read-on-a-subprocess-pipe-in-python

import sys
from subprocess import PIPE, Popen
from threading  import Thread

# magic numbers
# corresponds to enum values in pono/core/proverresult.h
UNKNOWN=255
FALSE=0
TRUE=1
ERROR=2

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run BMC and K-induction in parallel")
    parser.add_argument('btor_file')
    parser.add_argument('-k', '--bound', default='1000', help='The maximum bound to unroll to')
    parser.add_argument('-v', '--verbosity', action="store_true", help="Enable verbose output."
                        "   Note: this is buffered and only prints when a process finishes"
                        "         or there is an interrupt")

    args = parser.parse_args()
    btor_file = args.btor_file
    bound = args.bound
    verbosity = args.verbosity
    verbosity_option = '1' if args.verbosity else '0'
    pono="./pono"

    commands = {
        "BMC": [pono, '-e', 'bmc', '-v', verbosity_option, '-k', bound, btor_file],
        # "BMC+SimplePath": [pono, '-e', 'bmc-sp', '-v', verbosity_option, '-k', bound, btor_file],
        "K-Induction": [pono, '-e', 'ind', '-v', verbosity_option, '-k', bound, btor_file],
        "IC3": [pono, '-e', 'mbic3', '-v', verbosity_option, '-k', bound, btor_file],
        "ItpIC3": [pono, '-e', 'mbic3', '-v', verbosity_option, '-k', bound, '--ic3-indgen-mode', '2', btor_file],
        # give interpolant based methods a shorter bound -- impractical to go too large
        "Interpolant-based": [pono, '--smt-solver', 'msat', '-e', 'interp', '-v', verbosity_option, '-k', '100', btor_file],
        "ProphInterp-Arrays": [pono, '--smt-solver', 'msat', '-e', 'interp', '-v', verbosity_option, '-k', '100', '--ceg-prophecy-arrays', btor_file]
    }

    interp_processes = set()
    all_processes = []
    queues = {}
    name_map = {}

    # this one gets updated on the fly as processes end
    processes = []

    for name, cmd in commands.items():
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        if 'interp' in cmd:
            interp_processes.add(proc)
        processes.append(proc)
        all_processes.append(proc)
        name_map[proc] = name

    def print_process_output(proc):
        if proc.poll() is None:
            proc.terminate()
            proc.kill()
        out, _ = proc.communicate()
        print(out.decode('utf-8'))
        print()
        sys.stdout.flush()


    shutdown = False
    def handle_signal(signum, frame):
        # send signal recieved to subprocesses
        global shutdown
        if not shutdown:
            shutdown = True

            global verbosity
            if verbosity:
                for proc in all_processes:
                    print("{} output:".format(name_map[proc]))
                    print_process_output(proc)
                    print()
                    sys.stdout.flush()
            sys.exit(0)

    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)


    while not shutdown:
        for p in processes:
            if p.poll() is not None:
                # return code for unknown is 2
                # anything higher than that is an error
                # keep solving unless there are no more running processes
                if p.returncode >= 2:
                    processes.remove(p)
                    # print unknown only if this is the last process
                    if not processes:
                        print_process_output(p)
                        shutdown = True
                else:
                    # HACK don't return counter-examples from interpolation-based procedure
                    #      mathsat might return constant arrays in witness which can't be
                    #      printed in btor2 format
                    if p in interp_processes and p.returncode == FALSE:
                        processes.remove(p)
                        # this shouldn't happen but let's handle it just in case
                        if not processes:
                            print_process_output(bmc)
                            shutdown = True
                            break
                    else:
                        print_process_output(p)
                        shutdown = True
                        print("Winner was:", name_map[p])
                        break

        # just a double check
        if not processes:
            shutdown = True
            break

        if not shutdown:
            time.sleep(.001)

    # clean up
    for p in all_processes:
        if p.poll() is None:
            p.terminate()
