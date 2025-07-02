#include <utility>

#include "core/fts.h"
#include "core/prop.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "modifiers/mod_ts_prop.h"
#include "modifiers/prop_monitor.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

pair<Term, TransitionSystem> input_property_sys(SmtSolver & solver)
{
  Sort bvsort8 = solver->make_sort(BV, 8);
  FunctionalTransitionSystem fts(solver);
  Term max_val = fts.make_term(10, bvsort8);
  counter_system(fts, max_val);
  Term x = fts.named_terms().at("x");

  // add an input variable
  Term in = fts.make_inputvar("in", bvsort8);
  // constrain input to be less than a value
  fts.add_constraint(fts.make_term(BVUlt, in, fts.make_term(5, bvsort8)));

  Term prop_term = solver->make_term(
      BVUlt, solver->make_term(BVAdd, x, in), solver->make_term(15, bvsort8));
  return { prop_term, fts };
}

class PromoteInputvarsTests : public ::testing::Test,
                              public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override { s = create_solver(GetParam()); }
  SmtSolver s;
};

TEST_P(PromoteInputvarsTests, NoPromotion)
{
  auto res = input_property_sys(s);
  // need a property monitor
  TransitionSystem & ts = res.second;
  Term prop = add_prop_monitor(ts, res.first);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

TEST_P(PromoteInputvarsTests, WithPromotion)
{
  auto res = input_property_sys(s);
  // need a property monitor
  Term prop = res.first;
  TransitionSystem ts = promote_inputvars(res.second);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedPromoteInputvarsTests,
                         PromoteInputvarsTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
