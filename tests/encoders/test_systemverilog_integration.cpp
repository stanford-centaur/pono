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
  check_prover<KInduction>(
      "integration/arbiter_bus.sv", 12, ProverResult::TRUE);
}

// Integration #1b: the same design's liveness property, via the
// LivenessToSafetyTranslator path.
TEST_P(SVUnitTests, ArbiterBusLivenessFails)
{
  check_liveness_bmc("integration/arbiter_bus_live.sv", 8);
}

// Integration #2: register-file read-after-write invariant, proved
// via KInduction -- it is 1-inductive.
TEST_P(SVUnitTests, RegfileCtrlKInductionProves)
{
  check_prover<KInduction>(
      "integration/regfile_ctrl.sv", 6, ProverResult::TRUE);
}

// Integration #3: a credit-based link whose `assume property` keeps
// `push` from firing at zero credits. Without the assumption
// reaching the solver the counter wraps below zero and the bound is
// violated, so this is really a test that the assumption is there.
//
// Checked with IC3Bits, the one engine in this suite that is neither
// Bmc nor KInduction: it consumes the transition system differently
// and wants solver options of its own, so it is also a check that
// the frontend's output travels.
TEST_P(SVUnitTests, CreditLinkAssumeHonoredByIC3)
{
  check_prover<IC3Bits>("integration/credit_link.sv",
                        10,
                        ProverResult::TRUE,
                        /*filelists=*/{},
                        Engine::IC3_BITS);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVIntegrationTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
