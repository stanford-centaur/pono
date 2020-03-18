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

For more information, see the examples in `smt-switch/tests/btor`. Other useful files to visit include:
* `smt-switch/include/solver.h`: this is the main interface you will be using
* `smt-switch/include/ops.h`: this contains all the ops you might need
  * Note: create indexed ops like `Op(Extract, 7, 4)`

## Python bindings
To build the `cosa2` python bindings, first make sure that you have [Cython](https://cython.org/) version >= 0.29 installed. Then ensure that `smt-switch` and its python bindings are installed. Finally, you can configure with `./configure.sh --python` and then build normally. The sequence of commands would be as follows:
```
# Optional recommended step: start a python virtualenv
# If you install in the virtualenv, you will need to activate it each time before using pycosa2
# and deactivate the virtualenv with: deactivate
python3 -m venv env
source ./env/bin/activate
pip install Cython==0.29 pytest
./contrib/setup-smt-switch.sh --python
./contrib/setup-btor2tools.sh
pip install -e ./deps/smt-switch/build/python
./configure.sh --python
cd build
make -j4
pip install -e ./python
cd ../
# Test the bindings
pytest ./tests
```
