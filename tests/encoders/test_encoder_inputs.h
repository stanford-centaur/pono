#pragma once

// used in the individual tests to recover string paths from macros
#define STRHELPER(A) #A
#define STRFY(A) STRHELPER(A)

#include <cstddef>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "core/proverresult.h"
#include "smt-switch/solver_enums.h"

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

/** The query of a check-system command of a MoXI file, with the result of
 *  checking its property: FALSE if the query is satisfiable, TRUE if it is
 *  not, and UNKNOWN if k-induction cannot tell up to the bound the tests use.
 */
struct MoxiQuery
{
  string file;
  size_t index;  ///< of the check-system command
  unordered_set<smt::SolverAttribute> theories;  ///< what the solver needs
  pono::ProverResult result;
  unordered_set<smt::SolverEnum> excluded{};  ///< solvers that cannot tell
};

inline ostream & operator<<(ostream & os, const MoxiQuery & query)
{
  return os << query.file << "#" << query.index;
}

const vector<MoxiQuery> moxi_queries({
    { "assumptions.moxi", 0, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "assumptions.moxi", 1, { smt::THEORY_INT }, pono::ProverResult::TRUE },
    { "assumptions.moxi", 2, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "assumptions.moxi", 3, { smt::THEORY_INT }, pono::ProverResult::TRUE },
    { "constants.moxi", 0, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "constants.moxi", 1, { smt::THEORY_INT }, pono::ProverResult::UNKNOWN },
    { "constants.moxi", 2, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "counter.moxi", 0, { smt::THEORY_BV }, pono::ProverResult::FALSE },
    { "counter.moxi", 1, { smt::THEORY_BV }, pono::ProverResult::TRUE },
    { "counter.moxi", 2, { smt::THEORY_BV }, pono::ProverResult::FALSE },
    { "counter.moxi", 3, { smt::THEORY_BV }, pono::ProverResult::TRUE },
    { "current.moxi", 0, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "current.moxi", 1, { smt::THEORY_INT }, pono::ProverResult::TRUE },
    { "current.moxi", 2, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "enums.moxi",
      0,
      { smt::THEORY_INT, smt::THEORY_DATATYPE },
      pono::ProverResult::FALSE },
    { "enums.moxi",
      1,
      { smt::THEORY_INT, smt::THEORY_DATATYPE },
      pono::ProverResult::TRUE },
    { "enums.moxi",
      2,
      { smt::THEORY_INT, smt::THEORY_DATATYPE },
      pono::ProverResult::TRUE },
    { "enums.moxi",
      3,
      { smt::THEORY_INT, smt::THEORY_DATATYPE },
      pono::ProverResult::TRUE },
    // Bitwuzla answers unknown for the equality with the constant array.
    { "features.moxi",
      0,
      { smt::THEORY_BV, smt::CONSTARR },
      pono::ProverResult::FALSE,
      { smt::BZLA } },
    { "features.moxi",
      1,
      { smt::THEORY_BV, smt::CONSTARR },
      pono::ProverResult::TRUE,
      { smt::BZLA } },
    { "checks.moxi", 0, { smt::THEORY_BV }, pono::ProverResult::FALSE },
    { "checks.moxi", 1, { smt::THEORY_BV }, pono::ProverResult::TRUE },
    { "checks.moxi", 2, { smt::THEORY_BV }, pono::ProverResult::FALSE },
    { "reachables.moxi", 0, { smt::THEORY_BV }, pono::ProverResult::FALSE },
    { "reachables.moxi", 1, { smt::THEORY_BV }, pono::ProverResult::TRUE },
    { "reachables.moxi", 2, { smt::THEORY_BV }, pono::ProverResult::TRUE },
    // MathSAT makes the numeral 0 an integer even as a real, and smt-switch's
    // MathSAT backend refuses an ite with an integer and a real branch.
    { "reals.moxi",
      0,
      { smt::THEORY_REAL },
      pono::ProverResult::FALSE,
      { smt::MSAT } },
    { "reals.moxi",
      1,
      { smt::THEORY_REAL },
      pono::ProverResult::TRUE,
      { smt::MSAT } },
    { "real_arrays.moxi",
      0,
      { smt::THEORY_INT, smt::THEORY_REAL, smt::CONSTARR },
      pono::ProverResult::FALSE },
    { "real_arrays.moxi",
      1,
      { smt::THEORY_INT, smt::THEORY_REAL, smt::CONSTARR },
      pono::ProverResult::TRUE },
    { "repeated_args.moxi", 0, { smt::THEORY_INT }, pono::ProverResult::FALSE },
    { "repeated_args.moxi", 1, { smt::THEORY_INT }, pono::ProverResult::TRUE },
    { "repeated_args.moxi", 2, { smt::THEORY_INT }, pono::ProverResult::TRUE },
    { "subsystems.moxi", 0, { smt::THEORY_BV }, pono::ProverResult::FALSE },
    { "subsystems.moxi", 1, { smt::THEORY_BV }, pono::ProverResult::TRUE },
});

}  // namespace pono_tests
