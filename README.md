# cosa2
Next generation cosa.

## Setup

* Run `./contrib/setup-smt-switch.sh` -- it will build smt-switch with boolector and mathsat.
* Run `./contrib/setup-btor2tools.sh`.
* Run `./configure.sh`.
* Run `cd build`.
* Run `make`.

## Existing code

### Transition Systems
There are two Transition System interfaces:
* FunctionalTransitionSystem in fts.*
* RelationalTransitionSystem in rts.*


### Smt-Switch
Smt-switch is a C++ solver-agnostic API for SMT solvers. The main thing to remember is that everything is a pointer. Objects might be "typedef-ed" with `using` statements, but they're still `shared_ptr`s. Thus, when using a solver or a term, you need to use `->` accesses.

For more information, see the examples in `smt-switch/tests/btor`. Other useful files to visit include:
* `smt-switch/include/solver.h`: this is the main interface you will be using
* `smt-switch/include/ops.h`: this contains all the ops you might need
  * Note: create indexed ops like `Op(Extract, 7, 4)`
