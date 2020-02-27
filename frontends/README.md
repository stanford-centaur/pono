README

The main function is in `smv_cosa2.cpp` 

The parser needs latest version of bison(3.4.2 or newer) and flex(2.6.4) to be compatible with c++17

## Setup
Run `source set-env.sh` first then `make smvparser`
also need `otool -L ./smvparser`
`make clean` for cleaning all generated files
To parse a system run `./smvparser filename`
For example, run `./smvparser simple_counter.smv` 

## To do:
need to change some syntax parsing

## Not Include:
case expression
FAIRNESS
FUNCTION
MODULE
LTL
CTL 
PSL



