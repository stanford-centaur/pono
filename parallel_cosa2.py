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

    p0 = subprocess.Popen(['./cosa2', '-v', verbosity, '-k', bound, btor_file],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1 = subprocess.Popen(['./cosa2', '-v', verbosity, '-k', bound, '-i', btor_file],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    processes = [p0, p1]

    try:
        while processes:
            for i, p in enumerate(processes):
                if p.poll() is not None:
                    out, _ = p.communicate()
                    if out is not None:
                        print(out.decode('utf-8', errors='replace'))

                    # return code for unknown is 2 -- keep solving
                    if p.returncode == 2:
                        processes.remove(p)
                    else:
                        for pp in processes:
                            if pp != p:
                                pp.terminate()
                        processes = []
                        break

            if processes:
                time.sleep(.01)

    except KeyboardInterrupt as e:
        if args.verbosity:
            # on keyboard interrupt print the pipes if they haven't already been printed
            print("Got an interrupt, printing output")
            if p0.returncode is None:
                print("BMC output:")
                out, _ = p0.communicate()
                print(out.decode('utf-8', errors='replace'))

            if p1.returncode is None:
                print("K-induction output:")
                out, _ = p1.communicate()
                print(out.decode('utf-8', errors='replace'))
