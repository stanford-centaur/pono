#pragma once

// used in the individual tests to recover string paths from macros
#define STRHELPER(A) #A
#define STRFY(A) STRHELPER(A)

#include <string>
#include <unordered_map>
#include <vector>

#include "core/proverresult.h"

using namespace std;

namespace pono_tests {

// Named without the .btor2 extension that all of them share.
const vector<string> btor2_inputs({ "counter",
                                    "counter_true",
                                    "mem",
                                    "array_neq",
                                    "ridecore",
                                    "state2input",
                                    "write_counter" });

const vector<string> coreir_inputs({ "counters.json",
                                     "WrappedPE_nofloats.json",
                                     "SimpleALU.json" });

const unordered_map<string, pono::ProverResult> smv_inputs(
    { { "simple_counter.smv", pono::ProverResult::TRUE },
      { "simple_counter_integer.smv", pono::ProverResult::TRUE },
      { "simple_counter_integer_uf.smv", pono::ProverResult::FALSE },
      { "combined-false.smv", pono::ProverResult::FALSE },
      { "combined-true.smv", pono::ProverResult::TRUE },
      { "counter_bitvector.smv", pono::ProverResult::FALSE },
      { "counter_boolean.smv", pono::ProverResult::FALSE },
      { "signed_comparison.smv", pono::ProverResult::TRUE } });

}  // namespace pono_tests
