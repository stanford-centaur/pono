[![CI](https://github.com/stanford-centaur/pono/actions/workflows/ci.yml/badge.svg)](https://github.com/stanford-centaur/pono/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://github.com/stanford-centaur/pono/blob/main/LICENSE)

# Pono: A Flexible and Extensible SMT-Based Model Checker

Pono is a performant, adaptable, and extensible SMT-based model checker implemented in C++. It leverages [Smt-Switch](https://github.com/stanford-centaur/smt-switch), a generic C++ API for SMT solving. Pono was developed as the next
generation of [CoSA](https://github.com/cristian-mattarei/CoSA) and thus was originally named _cosa2_.

[Pono](http://wehewehe.org/gsdl2.85/cgi-bin/hdict?e=q-11000-00---off-0hdict--00-1----0-10-0---0---0direct-10-ED--4--textpukuielbert%2ctextmamaka-----0-1l--11-en-Zz-1---Zz-1-home-pono--00-3-1-00-0--4----0-0-11-10-0utfZz-8-00&a=d&d=D18537) is the Hawaiian word for proper, correct, or goodness. It is often used colloquially in the moral sense of "goodness" or "rightness," but also refers to "proper procedure" or "correctness." We use the word for multiple meanings. Our goal is that Pono can be a useful tool for people to verify the _correctness_ of systems, which is surely the _right_ thing to do.

## Publications

* Makai Mann, Ahmed Irfan, Florian Lonsing, Yahan Yang, Hongce Zhang, Kristopher Brown, Aarti Gupta, Clark W. Barrett: [Pono: A Flexible and Extensible SMT-Based Model Checker](https://link.springer.com/chapter/10.1007/978-3-030-81688-9_22). CAV 2021.
  * Evaluated [software artifact](https://figshare.com/articles/software/CAV_2021_Artifact_Pono_Model_Checker/14479542).
* Makai Mann, Amalee Wilson, Yoni Zohar, Lindsey Stuntz, Ahmed Irfan, Kristopher Brown, Caleb Donovick, Allison Guman, Cesare Tinelli, Clark W. Barrett: [Smt-Switch: A Solver-Agnostic C++ API for SMT Solving](https://link.springer.com/chapter/10.1007/978-3-030-80223-3_26). SAT 2021.
* Makai Mann: [Augmenting transition systems for scalable symbolic model checking](https://searchworks.stanford.edu/view/13972018). PhD thesis, Stanford University, 2021.

## Awards

Pono was awarded the Oski Award under its original name _cosa2_ at [HWMCC'19](http://fmv.jku.at/hwmcc19/) for solving the largest number of benchmarks overall.

## Setup

* [optional] Install bison and flex
  * If you don't have bison and flex installed globally, run `./contrib/setup-bison.sh` and `./contrib/setup-flex.sh`
  * Even if you do have bison, you might get errors about not being able to load `-ly`. In such a case, run the bison setup script.
* Run `./contrib/setup-smt-switch.sh` -- it will build smt-switch with Bitwuzla
  * [optional] to build with MathSAT (required for interpolation-based model checking) you need to obtain the libraries yourself
    * note that MathSAT is under a custom non-BSD compliant license and you must assume all responsibility for meeting the conditions
    * download the solver from https://mathsat.fbk.eu/download.html, unpack it and rename the directory to `./deps/mathsat`
    * then add the `--with-msat` flag to the `setup-smt-switch.sh` command.
* Run `./contrib/setup-btor2tools.sh`.
* Run `./configure.sh`.
  * if building with mathsat, also include `--with-msat` as an option to `configure.sh`
* Run `cd build`.
* Run `make`.
* [optional] Run `make check` to build and run the tests.

### Dependencies

* Please see the [README of smt-switch](https://github.com/stanford-centaur/smt-switch#dependencies) for required dependencies.
* Note to Arch Linux users: building Pono will fail if the static library of [GMP](https://gmplib.org/), which is required by [cvc5](https://github.com/cvc5/cvc5/blob/main/INSTALL.rst), is not installed on your system. You can fix it by installing `libgmp-static` from [AUR](https://aur.archlinux.org/packages/libgmp-static).

### Docker

We provide a Dockerfile for building a container image with Pono and all its dependencies.
To build the image locally, run:

```bash
docker build -t pono .
```

(You can replace `docker` with `podman` if you prefer a Docker-compatible alternative.)

Note that the Docker image does not include the MathSAT backend due to licensing restrictions.
If you need MathSAT support, please modify the Dockerfile manually to include it.

For your convenience, we publish the image built from the latest commit of main branch to GitHub Container Registry.
You can pull the image with:

```bash
docker pull ghcr.io/stanford-centaur/pono:latest
```

### Profiling

We link against the [gperftools library](https://github.com/gperftools/gperftools)
to generate profiling data. To enable profiling, run `./configure.sh`
with flag `--with-profiling` and recompile Pono by running `make` in
the `build` directory. This assumes that you have installed the
gperftools library before, e.g., on Ubuntu, run `sudo apt-get install
google-perftools libgoogle-perftools-dev`.

To profile, run `./pono --profiling-log=<log-file> ...` where
`<log-file>` is the path to a file where profiling data is
written. After normal or abnormal (e.g. via sending a signal)
termination of Pono, the collected profile data can be analyzed by
running, e.g., `google-pprof --text ./pono <log-file>` to produce a
textual profile. See `man google-pprof` for further options.

In general, see
[https://github.com/gperftools/gperftools/tree/master/docs](https://github.com/gperftools/gperftools/tree/master/docs)
for further information about gperftools.

gperftools is licensed under a BSD 3-clause license, see
[https://github.com/gperftools/gperftools/blob/master/COPYING](https://github.com/gperftools/gperftools/blob/master/COPYING).

## Existing code

### Transition Systems

There are two Transition System interfaces:

* FunctionalTransitionSystem in fts.*
* TransitionSystem in rts.*

### Smt-Switch

[Smt-switch](https://github.com/stanford-centaur/smt-switch) is a C++ solver-agnostic API for SMT solvers. The main thing to remember is that everything is a pointer. Objects might be "typedef-ed" with `using` statements, but they're still `shared_ptr`s. Thus, when using a solver or a term, you need to use `->` accesses.

For more information, see the example usage in the [smt-switch tests](https://github.com/stanford-centaur/smt-switch/tree/master/tests/btor).
Other useful files to visit include:

* `smt-switch/include/solver.h`: this is the main interface you will be using
* `smt-switch/include/ops.h`: this contains all the ops you might need
  * Note: create indexed ops like `Op(Extract, 7, 4)`

## Python bindings

To build the `pono` python bindings, first make sure that you have [Cython](https://cython.org/) version >= 3.0.0 installed. Then ensure that `smt-switch` and its python bindings are installed. Finally, you can configure with `./configure.sh --python` and then build normally. The bindings will be built under 'build/python', from where you can install them using `pip`. The sequence of commands would be as follows:

```bash
# Optional recommended step: start a python virtualenv
# If you install in the virtualenv, you will need to activate it each time before using pono
# and deactivate the virtualenv with: deactivate
python3 -m venv env
source ./env/bin/activate
# Install Python dependencies.
pip install Cython
# This will download and build smt-switch with bitwuzla and cvc5 in deps/smt-switch.
./contrib/setup-smt-switch.sh --python
# Install the smt-switch Python bindings in the virtualenv.
pip install ./deps/smt-switch/build/python
# Download and build the btor2tools dependency.
./contrib/setup-btor2tools.sh
# Configure and build Pono.
./configure.sh --python
make -C build -j4  # change the `-j4` to match the number of CPU cores on your system.
# Install the built Python bindings in the virtualenv (including the `pytest` dependency for testing).
pip install './build/python[test]'
# Test the bindings.
pytest ./tests
```

## Documentation

To generate documentation from the C++ source files, install [Doxygen](https://www.doxygen.nl/index.html), configure with `./configure.sh --docs`, then build the `docs` target:

```
# Install Doxygen
# Ubuntu
sudo apt-get install doxygen
# Arch Linux
sudo pacman -S doxygen
# macOS
brew install doxygen

# Build documentation
./configure.sh --docs
cmake --build build -t docs

# Open the HTML index
# Linux
xdg-open build/html/index.html
# macOS
open build/html/index.html
```

## Generating BTOR2 from Verilog

The best tool for creating BTOR2 from Verilog is [Yosys](https://github.com/YosysHQ/yosys). Yosys has an excellent manual [here](http://www.clifford.at/yosys/files/yosys_manual.pdf). 
You can also run yosys interactively by running yosys with no arguments. Then you can view help messages for each command with: `help <command>`. Running `help` with no arguments lists all commands.

A particularly useful command if you're having trouble is `show`, which can show the current state of the circuit in Yosys.

### Yosys Quick Start

Below is an example file with comments explaining each command that produces a BTOR2 file for [./samples/counter-false.v](./samples/counter-false.v). This should be enough for most use cases.

Once you have `yosys` installed, copy the text below into `gen-btor.ys` in the top-level of this repository. Then, running `yosys -s gen-btor.ys` will produce the BTOR2 file.

```
# read in the file(s) -- there can be multiple
# whitespace separated files, and you can
# escape new lines if necessary
read -formal ./samples/counter-false.v;

# prep does a conservative elaboration
# of the top module provided
prep -top counter;

# this command just does a sanity check
# of the hierarchy
hierarchy -check;

# If an assumption is flopped, you might
# see strange behavior at the last state
# (because the clock hasn't toggled)
# this command ensures that assumptions
# hold at every state
chformal -assume -early;

# this processes memories
# nomap means it will keep them as arrays
memory -nomap;

# flatten the design hierarchy
flatten;

# (optional) uncomment and set values to simulate reset signal
# use -resetn for an active low pin
# -n configures the number of cycles to simulate
# -rstlen configures how long the reset is active (recommended to keep it active for the whole simulation)
# -w tells it to write back the final state of the simulation as the initial state in the btor2 file
# another useful option is -zinit which zero initializes any uninitialized state
# sim -clock <clockpin> -reset <resetpin> -n <number of cycles> -rstlen <number of cycles> -w <top_module>

# (optional) use an "explicit" clock
# e.g. every state is a half cycle of the
# fastest clock
# use this option if you see errors that
# refer to "adff" or asynchronous components
# IMPORTANT NOTE: the clocks are not
# automatically toggled if you use this option
# clk2fflogic;

# This turns all undriven signals into
# inputs
setundef -undriven -expose;

# This writes to a file in BTOR2 format
write_btor counter-false.btor2
```
