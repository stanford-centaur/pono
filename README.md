# cosa2
Next generation cosa.

## Setup

* Run `./contrib/setup-smt-switch.sh` -- it will build smt-switch with boolector
  * add `--with-msat` to also build with MathSAT (required for interpolant-based model checking)
  * note that MathSAT is under a custom non-BSD compliant license and you must assume all responsibility for meeting the conditions
* Run `./contrib/setup-btor2tools.sh`.
* Run `./configure.sh`.
* Run `cd build`.
* Run `make`.

## Existing code

### Transition Systems
There are two Transition System interfaces:
* FunctionalTransitionSystem in fts.*
* TransitionSystem in rts.*


### Smt-Switch
Smt-switch is a C++ solver-agnostic API for SMT solvers. The main thing to remember is that everything is a pointer. Objects might be "typedef-ed" with `using` statements, but they're still `shared_ptr`s. Thus, when using a solver or a term, you need to use `->` accesses.

For more information, see the example usage in the [smt-switch tests](https://github.com/makaimann/smt-switch/tree/master/tests/btor).
Other useful files to visit include:
* `smt-switch/include/solver.h`: this is the main interface you will be using
* `smt-switch/include/ops.h`: this contains all the ops you might need
  * Note: create indexed ops like `Op(Extract, 7, 4)`

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

# If an assumption is flopped, you might
# see strange behavior at the last state
# (because the clock hasn't toggled)
# this command ensures that assumptions
# hold at every state
chformal -assume -early;

# This turns all undriven signals into
# inputs
setundef -undriven -expose;

# This writes to a file in BTOR2 format
write_btor counter-false.btor2
```
