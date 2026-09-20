// Tests that drive the pono executable itself, for justice properties and the
// fairness constraints that go with them.

#include <string>

#include "cli_test_fixture.h"
#include "gtest/gtest.h"

using namespace std;

namespace pono_tests {

class JusticeCliUnitTests : public CliUnitTests
{
};

// A justice property with no fairness constraint: the lasso at mode=1
// satisfies the condition infinitely often, so it is a counterexample.
TEST_F(JusticeCliUnitTests, JusticeOnly)
{
  const PonoRun run =
      run_pono({ "--justice", input_path("btor2/justice_only.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// The same property under the other translator, which counts observations of
// the single condition rather than translating it to a safety property.
TEST_F(JusticeCliUnitTests, JusticeOnlyWithKLiveness)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--justice-translator",
                                 "klive",
                                 input_path("btor2/justice_only.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// A justice property that holds, which is the outcome the cases above cannot
// reach: no lasso satisfies the condition infinitely often.
TEST_F(JusticeCliUnitTests, JusticeHolds)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--engine",
                                 "ind",
                                 "--bound",
                                 "20",
                                 input_path("btor2/justice_holds.btor2") });
  EXPECT_EQ(run.output, "unsat\nj0\n");
}

// Each justice line is a separate property, so --prop picks between them. The
// two here disagree, which is what makes the selection visible.
TEST_F(JusticeCliUnitTests, SecondJusticeProperty)
{
  const PonoRun run =
      run_pono({ "--justice",
                 "--prop",
                 "1",
                 "--engine",
                 "ind",
                 "--bound",
                 "20",
                 input_path("btor2/justice_two_properties.btor2") });
  EXPECT_EQ(run.output, "unsat\nj1\n");
}

// Reducing to the state the justice condition depends on leaves the verdict
// alone.
TEST_F(JusticeCliUnitTests, JusticeWithConeOfInfluence)
{
  const PonoRun run =
      run_pono({ "--justice",
                 "--static-coi",
                 input_path("btor2/justice_unrelated_state.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// k-liveness takes its own cone of influence rather than going through the
// one the driver takes, so the reduction reaches it by a separate route.
TEST_F(JusticeCliUnitTests, JusticeWithConeOfInfluenceAndKLiveness)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--justice-translator",
                                 "klive",
                                 "--static-coi",
                                 input_path("btor2/justice_only.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// The justice condition alone is violated, so the property only holds if the
// fair line is honored. k-induction is needed because bmc, the default
// engine, cannot prove a property.
TEST_F(JusticeCliUnitTests, JusticeHonorsFairnessConstraints)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--engine",
                                 "ind",
                                 input_path("btor2/fair_frozen_mode.btor2") });
  EXPECT_TRUE(contains(run.output, "unsat\nj0"));
}

// Btor2 fairness constraints apply to justice properties only, so a fair line
// has no bearing on checking a bad property. The whole output is matched
// rather than searched, because "sat" occurs in "unsat" too.
TEST_F(JusticeCliUnitTests, FairnessDoesNotAffectSafetyChecking)
{
  const PonoRun run =
      run_pono({ input_path("btor2/fair_with_bad_property.btor2") });
  EXPECT_EQ(run.output, "sat\nb0\n");
}

// The cone of influence is taken over the fairness constraints as well, so
// the state only they watch survives it. Reducing over the justice condition
// alone drops that state and the run dies translating the property.
TEST_F(JusticeCliUnitTests, JusticeWithConeOfInfluenceKeepsFairnessState)
{
  const PonoRun run =
      run_pono({ "--justice",
                 "--static-coi",
                 input_path("btor2/fair_independent_state.btor2") });
  EXPECT_EQ(run.output, "sat\nj0\n");
}

// k-liveness counts one condition, so it has to reject a fairness constraint
// rather than silently check the justice condition alone.
TEST_F(JusticeCliUnitTests, KLivenessRejectsFairness)
{
  // Matched on enough of the message to tell it apart from a path, which can
  // itself contain the word "fairness".
  expect_rejected(run_pono({ "--justice",
                             "--justice-translator",
                             "klive",
                             input_path("btor2/fair_frozen_mode.btor2") }),
                  "including fairness constraints");
}

// Errors label the property by the kind that was requested, so a failure
// under --justice reports j, not b.
TEST_F(JusticeCliUnitTests, ErrorLabelsPropertyByRequestedKind)
{
  const PonoRun run = run_pono({ "--justice",
                                 "--prop",
                                 "5",
                                 input_path("btor2/fair_frozen_mode.btor2") });
  expect_rejected(run, "Property index 5");
#ifdef NDEBUG
  // Only the build that catches the exception reaches the result printing.
  EXPECT_TRUE(contains(run.output, "error\nj5"));
#endif
}

// The SMV and VMT encoders parse no liveness properties, so --justice has to
// be refused instead of silently checking an invariant property.
TEST_F(JusticeCliUnitTests, JusticeRejectedForSmv)
{
  expect_rejected(run_pono({ "--justice", input_path("smv/counter.smv") }),
                  "--justice is not supported for smv");
}

}  // namespace pono_tests
