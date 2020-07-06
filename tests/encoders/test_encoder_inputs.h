#pragma once

// used in the individual tests to recover string paths from macros
#define STRHELPER(A) #A
#define STRFY(A) STRHELPER(A)

#include <string>
#include <vector>

using namespace std;

namespace pono_tests {

const vector<string> coreir_inputs({ "counters.json",
                                     "WrappedPE_nofloats.json",
                                     "SimpleALU.json" });

const vector<string> btor2_inputs({ "counter.btor",
                                    "counter-true.btor",
                                    "mem.btor",
                                    "ridecore.btor",
                                    "state2input.btor" });

}  // namespace pono_tests
