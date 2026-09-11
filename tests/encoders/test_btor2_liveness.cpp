#include <string>

#include "core/fts.h"
#include "core/prop.h"
#include "core/proverresult.h"
#include "engines/bmc.h"
#include "engines/kinduction.h"
#include "engines/kliveness.h"
#include "frontends/btor2_encoder.h"
#include "gtest/gtest.h"
#include "modifiers/liveness_to_safety_translator.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class Btor2LivenessUnitTests : public ::testing::Test,
                               public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string input_path(const string & name) const
  {
    return string(STRFY(PONO_SRC_DIR)) + "/tests/encoders/inputs/btor2/" + name;
  }

  /** The condition set a generalized-Buchi search is given for the property
   *  at index 0: its justice conditions plus the file's fairness constraints,
   *  matching how pono.cpp combines them.
   */
  TermVec all_conditions(const BTOR2Encoder & be) const
  {
    TermVec conditions = be.justicevec().at(0);
    conditions.insert(
        conditions.end(), be.fairvec().begin(), be.fairvec().end());
    return conditions;
  }
};

// Btor2 fair operands are always bitvectors of width 1, but the liveness to
// safety translation combines conditions with boolean operators, so the
// encoder has to convert them the same way it converts justice conditions.
TEST_P(Btor2LivenessUnitTests, FairTermsAreBoolSorted)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("fair_frozen_mode.btor2"), fts);
  ASSERT_EQ(be.fairvec().size(), 1);
  EXPECT_EQ(be.fairvec()[0]->get_sort()->get_sort_kind(), SortKind::BOOL);
}

// The justice condition alone is violated by the lasso at mode=1, but that
// lasso never satisfies the fairness constraint, so adding the constraint
// removes the only counterexample. This is what makes the union of the two
// condition sets more than a no-op.
TEST_P(Btor2LivenessUnitTests, FairnessRestrictsCounterexample)
{
  const string filename = input_path("fair_frozen_mode.btor2");
  const int bound = 10;

  SmtSolver justice_solver = create_solver(GetParam());
  justice_solver->set_opt("incremental", "true");
  FunctionalTransitionSystem justice_fts(justice_solver);
  BTOR2Encoder justice_be(filename, justice_fts);
  Term justice_term = LivenessToSafetyTranslator{}.translate(
      justice_fts, justice_be.justicevec().at(0));
  SafetyProperty justice_prop(justice_solver, justice_term);
  Bmc justice_bmc(justice_prop, justice_fts, justice_solver);
  EXPECT_EQ(justice_bmc.check_until(bound), ProverResult::FALSE);

  // A separate solver, because the translation adds equally named variables.
  SmtSolver fair_solver = create_solver(GetParam());
  fair_solver->set_opt("incremental", "true");
  FunctionalTransitionSystem fair_fts(fair_solver);
  BTOR2Encoder fair_be(filename, fair_fts);
  Term fair_term =
      LivenessToSafetyTranslator{}.translate(fair_fts, all_conditions(fair_be));
  SafetyProperty fair_prop(fair_solver, fair_term);
  Bmc fair_bmc(fair_prop, fair_fts, fair_solver);
  EXPECT_NE(fair_bmc.check_until(bound), ProverResult::FALSE);

  // And no deeper counterexample exists either.
  SmtSolver ind_solver = create_solver(GetParam());
  ind_solver->set_opt("incremental", "true");
  FunctionalTransitionSystem ind_fts(ind_solver);
  BTOR2Encoder ind_be(filename, ind_fts);
  Term ind_term =
      LivenessToSafetyTranslator{}.translate(ind_fts, all_conditions(ind_be));
  SafetyProperty ind_prop(ind_solver, ind_term);
  KInduction ind(ind_prop, ind_fts, ind_solver);
  EXPECT_EQ(ind.check_until(bound), ProverResult::TRUE);
}

// Fairness constraints that cannot all recur exclude every trace, so the
// property holds for want of any fair counterexample. Pono reports this as a
// plain proof, the same way it treats an unsatisfiable set of constraint
// lines in safety checking.
TEST_P(Btor2LivenessUnitTests, ContradictoryFairnessProvesVacuously)
{
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("fair_contradictory.btor2"), fts);
  ASSERT_EQ(be.fairvec().size(), 2);
  Term prop_term =
      LivenessToSafetyTranslator{}.translate(fts, all_conditions(be));
  SafetyProperty prop(s, prop_term);
  KInduction ind(prop, fts, s);
  EXPECT_EQ(ind.check_until(10), ProverResult::TRUE);
}

// k-liveness counts observations of a single condition, so it cannot honor a
// fairness constraint alongside a justice condition. It has to say so rather
// than check the justice condition on its own.
TEST_P(Btor2LivenessUnitTests, KLivenessRejectsFairness)
{
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("fair_frozen_mode.btor2"), fts);
  ASSERT_EQ(be.justicevec().at(0).size(), 1);
  ASSERT_EQ(be.fairvec().size(), 1);
  LivenessProperty prop(s, all_conditions(be));
  try {
    KLiveness kliveness(prop, fts, s, PonoOptions());
    FAIL() << "expected KLiveness to reject the fairness constraint";
  }
  catch (const PonoException & e) {
    EXPECT_NE(string(e.what()).find("fairness"), string::npos)
        << "message should name fairness constraints, got: " << e.what();
  }
}

// Without a fair line the condition set is just the justice conditions, so
// the counterexample the fairness constraint ruled out above is found again.
TEST_P(Btor2LivenessUnitTests, JusticeOnlyIsUnaffected)
{
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("justice_only.btor2"), fts);
  EXPECT_TRUE(be.fairvec().empty());
  EXPECT_EQ(all_conditions(be), be.justicevec().at(0));
  Term prop_term =
      LivenessToSafetyTranslator{}.translate(fts, all_conditions(be));
  SafetyProperty prop(s, prop_term);
  Bmc bmc(prop, fts, s);
  EXPECT_EQ(bmc.check_until(10), ProverResult::FALSE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverBtor2LivenessUnitTests,
                         Btor2LivenessUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
