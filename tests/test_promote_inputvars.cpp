#include <utility>

#include "core/fts.h"
#include "core/prop.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "modifiers/mod_ts_prop.h"
#include "modifiers/prop_monitor.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"

using namespace pono;
using namespace smt;
using namespace std;

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
  ts = promote_inputvars(ts, ivs_in_prop);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

TEST_P(PromoteInputvarsTests, PromoteAllInputs)
{
  ts = promote_inputvars(ts);

  SafetyProperty p(s, prop);
  KInduction kind(p, ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedPromoteInputvarsTests,
    PromoteInputvarsTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(std::vector<TSEnum>{ Functional,
                                                            Relational })));

}  // namespace pono_tests
