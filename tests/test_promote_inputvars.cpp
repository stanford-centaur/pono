#include <tuple>
#include <vector>

#include "core/fts.h"
#include "core/prop.h"
#include "core/proverresult.h"
#include "core/rts.h"
#include "core/ts.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "modifiers/mod_ts_prop.h"
#include "modifiers/prop_monitor.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

enum TSEnum
{
  Functional,
  Relational
};

Term input_property_sys(TransitionSystem & ts)
{
  Sort bvsort8 = ts.make_sort(BV, 8);
  Term max_val = ts.make_term(10, bvsort8);
  counter_system(ts, max_val);
  Term x = ts.named_terms().at("x");

  // add an input variable
  Term in = ts.make_inputvar("in", bvsort8);
  // constrain input to be less than a value
  ts.add_constraint(ts.make_term(BVUlt, in, ts.make_term(5, bvsort8)));

  Term prop_term = ts.make_term(
      BVUlt, ts.make_term(BVAdd, x, in), ts.make_term(15, bvsort8));
  return prop_term;
}

class PromoteInputvarsTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, TSEnum>>
{
 protected:
  void SetUp() override
  {
    SolverEnum se = std::get<0>(GetParam());
    s = create_solver(se, se == BTOR);
    if (std::get<1>(GetParam()) == Functional) {
      ts = FunctionalTransitionSystem(s);
    } else {
      ts = RelationalTransitionSystem(s);
    }
    prop = input_property_sys(ts);
  }
  SmtSolver s;
  Term prop;
  TransitionSystem ts;
};

TEST_P(PromoteInputvarsTests, AddPropMonitor)
{
  // need a property monitor
  prop = add_prop_monitor(ts, prop);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

TEST_P(PromoteInputvarsTests, PromoteInputsInProp)
{
  UnorderedTermSet ivs_in_prop;
  get_free_symbolic_consts(prop, ivs_in_prop);
  // the monitor has to come first: a property reading an input holds on the
  // transitions that input labels, so promoting it and then checking the
  // promoted property in every state would check one state too many
  prop = add_prop_monitor(ts, prop);
  ts = promote_inputvars(ts, ivs_in_prop);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

TEST_P(PromoteInputvarsTests, PromoteAllInputs)
{
  prop = add_prop_monitor(ts, prop);
  ts = promote_inputvars(ts);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

TEST_P(PromoteInputvarsTests, PromotionKeepsConstraintWindow)
{
  // the constraint reads an input, so it is enforced on the transitions that
  // input labels and not in the initial state
  ASSERT_EQ(ts.constraints().size(), 1);
  EXPECT_FALSE(ts.constraints().at(0).second);
  Term init_before = ts.init();

  // promoting the input turns it into a state variable, but that must not
  // widen the constraint to hold in the initial state as well
  TransitionSystem promoted = promote_inputvars(ts);
  EXPECT_TRUE(promoted.inputvars().empty());
  ASSERT_EQ(promoted.constraints().size(), 1);
  EXPECT_FALSE(promoted.constraints().at(0).second);
  EXPECT_EQ(promoted.init(), init_before);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedPromoteInputvarsTests,
    PromoteInputvarsTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(std::vector<TSEnum>{ Functional,
                                                            Relational })));

}  // namespace pono_tests
