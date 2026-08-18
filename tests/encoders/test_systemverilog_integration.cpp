#ifdef WITH_SLANG

#include "engines/ic3bits.h"
#include "engines/kinduction.h"
#include "options/options.h"
#include "sv_test_fixture.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

// ---------------------------------------------------------------------------
// Three self-contained multi-module designs, each combining several
// features already covered in isolation elsewhere in this suite, each
// checked with a *different* prover engine (every other test in this
// suite -- unit or otherwise -- only ever uses Bmc).
// ---------------------------------------------------------------------------

// Integration #1a: round-robin arbiter safety property, plain Bmc.
TEST_P(SVUnitTests, ArbiterBusSafetyHolds)
{
  check_bmc("integration/arbiter_bus.sv", 6, ProverResult::UNKNOWN);
}

// Integration #1b: the same design's liveness property, via the
// LivenessToSafetyTranslator path.
TEST_P(SVUnitTests, ArbiterBusLivenessFails)
{
  check_liveness_bmc("integration/arbiter_bus_live.sv", 8);
}

// Integration #2a: register-file read-after-write invariant, plain
// Bmc -- can only ever report UNKNOWN, not an inductive proof.
TEST_P(SVUnitTests, RegfileCtrlBmcInconclusive)
{
  check_bmc("integration/regfile_ctrl.sv", 6, ProverResult::UNKNOWN);
}

// Integration #2b: the same invariant, actually *proved* via
// KInduction (it is 1-inductive) -- the first test in this suite to
// run a non-Bmc engine over the SV frontend's output.
TEST_P(SVUnitTests, RegfileCtrlKInductionProves)
{
  check_prover<KInduction>(
      "integration/regfile_ctrl.sv", 6, ProverResult::TRUE);
}

// Integration #3: credit-based link, checked with IC3Bits.  Written
// by hand rather than via check_prover() because IC3's proof isn't
// bounded by a `check_until` step count the way Bmc's counterexample
// search is.  FIXED: this used to demonstrate the assume_property.sv
// gap (test_systemverilog_sva.cpp) mattering at design scale -- IC3
// found credits overflowing past 4'd4 by driving `push` while
// credits == 0, exactly the case `assume property` was meant to rule
// out but (before the ConcurrentAssertion fix) silently didn't.  Now
// that assumptions are honored, IC3 actually *proves* the property.
TEST_P(SVUnitTests, CreditLinkAssumeHonoredByIC3)
{
  // IC3Bits needs solver options (e.g. unsat-assumptions production)
  // beyond the plain incremental/produce-models pair Bmc/KInduction
  // get by away with -- create_solver_for(..., Engine::IC3_BITS, ...)
  // is the same helper test_ic3bits.cpp uses for exactly this reason.
  SmtSolver s = create_solver_for(GetParam(), Engine::IC3_BITS, false);
  FunctionalTransitionSystem fts(s);
  SystemVerilogEncoder enc(sv_path("integration/credit_link.sv"), fts);
  ASSERT_EQ(enc.propvec().size(), 1u);
  Term prop_term = enc.propvec()[0];

  TransitionSystem ts = fts;
  if (Term rst = find_reset(ts)) {
    Term reset_done = add_reset_seq(ts, rst, /*reset_bnd=*/1);
    prop_term = ts.solver()->make_term(Implies, reset_done, prop_term);
  }

  SafetyProperty prop(ts.solver(), prop_term);
  IC3Bits ic3bits(prop, ts, s);
  EXPECT_EQ(ic3bits.check_until(10), ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVIntegrationTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
