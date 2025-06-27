#include "core/fts.h"
#include "core/prop.h"
#include "core/rts.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "modifiers/mod_ts_prop.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class PseudoInitPropUnitTests : public ::testing::Test,
                                public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam(), false);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort bvsort8;
};

TEST_P(PseudoInitPropUnitTests, CounterSystemSafe)
{
  FunctionalTransitionSystem fts(s);
  counter_system(fts, fts.make_term(10, bvsort8));
  Term x = fts.named_terms().at("x");

  Term prop_term = s->make_term(BVUle, x, s->make_term(10, bvsort8));

  TransitionSystem rts = pseudo_init_and_prop(fts, prop_term);
  assert(!rts.is_functional());
  SafetyProperty p(s, prop_term);

  KInduction kind(p, rts, s);
  ProverResult r = kind.check_until(11);
  ASSERT_EQ(r, TRUE);
}

TEST_P(PseudoInitPropUnitTests, TrivialUnsafe)
{
  // need to be careful when modifying transition system and property
  // this trivial unsafe system is actually an edge case
  // that can be hard to be sure we avoid with this transformation

  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort8);
  // no initial state so property is trivially false
  Term prop_term = s->make_term(BVUle, x, s->make_term(10, bvsort8));

  // but has a transition relation that is empty (deadlocked)
  rts.set_trans(rts.make_term(false));

  TransitionSystem ts = pseudo_init_and_prop(rts, prop_term);
  assert(!ts.is_functional());
  SafetyProperty p(s, prop_term);

  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(3);
  ASSERT_EQ(r, FALSE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverPseudoInitPropUnitTests,
                         PseudoInitPropUnitTests,
                         testing::ValuesIn(available_solver_enums()));
}  // namespace pono_tests
