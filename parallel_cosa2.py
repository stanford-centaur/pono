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

# try:
#     from queue import Queue, Empty
# except ImportError:
#     from Queue import Queue, Empty  # python 2.x

ON_POSIX = 'posix' in sys.builtin_module_names

def enqueue_output(out, l):
    for line in iter(out.readline, b''):
        l.append(line.decode('utf-8'))
    # out.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run BMC and K-induction in parallel")
    parser.add_argument('btor_file')
    parser.add_argument('-k', '--bound', default='200', help='The maximum bound to unroll to')
    parser.add_argument('-v', '--verbosity', action="store_true", help="Enable verbose output."
                        "   Note: this is buffered and only prints when a process finishes"
                        "         or there is an interrupt")

    args = parser.parse_args()
    btor_file = args.btor_file
    bound = args.bound
    verbosity = args.verbosity
    verbosity_option = '1' if args.verbosity else '0'

    commands = {
        "BMC": ['./cosa2', '-e', 'bmc', '-v', verbosity_option, '-k', bound, btor_file],
        "BMC+SimplePath": ['./cosa2', '-e', 'bmc-sp', '-v', verbosity_option, '-k', bound, btor_file],
        "K-Induction": ['./cosa2', '-e', 'ind', '-v', verbosity_option, '-k', bound, btor_file],
        "Interpolant-based": ['./cosa2', '-e', 'interp', '-v', verbosity_option, '-k', bound, btor_file]
    }

    interp_processes = set()
    all_processes = []
    queues = {}
    name_map = {}

    # this one gets updated on the fly as processes end
    processes = []

    for name, cmd in commands.items():
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, close_fds=ON_POSIX)
        if 'interp' in cmd:
            interp_processes.add(proc)
        processes.append(proc)
        all_processes.append(proc)
        name_map[proc] = name

        # instantiate watcher thread
        q = []
        t = Thread(target=enqueue_output, args=(proc.stdout, q))
        queues[proc] = q
        t.daemon = True
        t.start()


    def print_process_output(proc):
        for line in queues[proc]:
            print(line, end='')
        print()
        sys.stdout.flush()


    shutdown = False
    def handle_signal(signum, frame):
        # send signal recieved to subprocesses
        global shutdown
        if not shutdown:
            shutdown = True
            for proc in processes:
                if proc.poll() is None:
                    proc.terminate()

            global verbosity
            if verbosity:
                # too slow to communicate with process in signal handling
                # use cached lines
                for proc in all_processes:
                    print("{} output:".format(name_map[proc]))
                    for line in queues[proc]:
                        print(line, end='')
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
                    if p in interp_processes and p.returncode == 1:
                        processes.remove(p)
                        # this shouldn't happen but let's handle it just in case
                        if not processes:
                            print_process_output(bmc)
                            shutdown = True
                            break
                    else:
                        print_process_output(p)
                        shutdown = True
                        break

        # just a double check
        if not processes:
            shutdown = True
            break

        if not shutdown:
            time.sleep(.001)
