#!/usr/bin/env python

import argparse
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
    verbosity = '1' if args.verbosity else '0'

    bmc = subprocess.Popen(['./cosa2', '-e', 'bmc', '-v', verbosity, '-k', bound, btor_file],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    induc = subprocess.Popen(['./cosa2', '-e', 'ind', '-v', verbosity, '-k', bound, btor_file],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    interp = subprocess.Popen(['./cosa2', '-e', 'interp', '-v', verbosity, '-k', bound, btor_file],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    interp_processes = set({interp})
    processes = [bmc, induc, interp]

    try:
        while processes:
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
                    else:
                        # HACK don't return counter-examples from interpolation-based procedure
                        #      mathsat might return constant arrays in witness which can't be
                        #      printed in btor2 format
                        if p in interp_processes and p.returncode == 1:
                            processes.remove(p)
                            # this shouldln't happen but let's handle it just in case
                            if not processes:
                                out, _ = bmc.communicate()
                                print(out.decode('utf-8', errors='replace'))
                                break
                        else:
                            out, _ = p.communicate()
                            if out is not None:
                                print(out.decode('utf-8', errors='replace'))
                            for pp in processes:
                                if pp != p:
                                    pp.terminate()
                            processes = []
                            break

            if processes:
                time.sleep(.01)

    except KeyboardInterrupt as e:
        if args.verbosity:
            # on keyboard interrupt, print the pipes
            # nothing should have printed yet, otherwise we would have already exited
            print("Got an interrupt, printing output")

            print("BMC output:")
            out, _ = bmc.communicate()
            print(out.decode('utf-8', errors='replace'))

            print("K-induction output:")
            out, _ = induc.communicate()
            print(out.decode('utf-8', errors='replace'))

            print("Interpolation-based model checking output:")
            out, _ = interp.communicate()
            print(out.decode('utf-8', errors='replace'))
