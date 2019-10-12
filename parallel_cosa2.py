#!/usr/bin/env python3

import argparse
import signal
import subprocess
import sys
import time

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

    bmc = subprocess.Popen(['./cosa2', '-e', 'bmc', '-v', verbosity_option, '-k', bound, btor_file],
                           stdout=subprocess.PIPE)
    bmc_sp = subprocess.Popen(['./cosa2', '-e', 'bmc-sp', '-v', verbosity_option, '-k', bound, btor_file],
                              stdout=subprocess.PIPE)
    induc = subprocess.Popen(['./cosa2', '-e', 'ind', '-v', verbosity_option, '-k', bound, btor_file],
                          stdout=subprocess.PIPE)
    interp = subprocess.Popen(['./cosa2', '-e', 'interp', '-v', verbosity_option, '-k', bound, btor_file],
                          stdout=subprocess.PIPE)

    shutdown = False

    # these stay constant
    interp_processes = set({interp})
    all_processes = [bmc, bmc_sp, induc, interp]

    # this one is updated
    processes = [bmc, bmc_sp, induc, interp]
    pos_map = {0: "BMC:",
               1: "BMC+SimplePath",
               2: "K-Induction",
               3: "Interpolant-based"}

    def handle_signal(signum, frame):
        # send signal recieved to subprocesses
        global shutdown
        if not shutdown:
            shutdown = True
            for i, proc in enumerate(all_processes):
                if proc.poll() is None:
                    proc.terminate()

                global verbosity
                if verbosity:
                    # this is too slow
                    # out, err = proc.communicate()
                    # print(out.decode('utf-8', errors='replace'))

                    print("{} output:".format(pos_map[i]))
                    for line in proc.stdout:
                        print(line.decode('utf-8'), end='')
                    print()
            sys.stdout.flush()
            sys.exit(0)

    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)


    while not shutdown:
        for i, p in enumerate(processes):
            if p.poll() is not None:
                # return code for unknown is 2
                # anything higher than that is an error
                # keep solving unless there are no more running processes
                if p.returncode >= 2:
                    processes.remove(p)
                    # print unknown only if this is the last process
                    if not processes:
                        out, _ = p.communicate()
                        if out is not None:
                            print(out.decode('utf-8', errors='replace'))
                            sys.stdout.flush()
                        shutdown = True
                else:
                    # HACK don't return counter-examples from interpolation-based procedure
                    #      mathsat might return constant arrays in witness which can't be
                    #      printed in btor2 format
                    if p in interp_processes and p.returncode == 1:
                        processes.remove(p)
                        # this shouldn't happen but let's handle it just in case
                        if not processes:
                            out, _ = bmc.communicate()
                            print(out.decode('utf-8', errors='replace'))
                            sys.stdout.flush()
                            shutdown = True
                            break
                    else:
                        out, _ = p.communicate()
                        if out is not None:
                            print(out.decode('utf-8', errors='replace'))
                            sys.stdout.flush()
                        for pp in processes:
                            if pp != p:
                                pp.terminate()
                        processes = []
                        shutdown = True
                        break

        # just a double check
        if not processes:
            shutdown = True

        if not shutdown:
            time.sleep(.01)
