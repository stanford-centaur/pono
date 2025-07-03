#include "core/fts.h"
#include "core/functional_unroller.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "gtest/gtest.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class UnrollerUnitTests : public ::testing::Test,
                          public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    bvsort = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort bvsort;
};

TEST_P(UnrollerUnitTests, FTS_Unroll)
{
  FunctionalTransitionSystem fts(s);
  counter_system(fts, fts.make_term(10, bvsort));
  Term x = fts.named_terms().at("x");

  Unroller u(fts);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
}

TEST_P(UnrollerUnitTests, RTS_Unroll)
{
  RelationalTransitionSystem rts(s);
  counter_system(rts, rts.make_term(10, bvsort));
  Term x = rts.named_terms().at("x");

  Unroller u(rts);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
  Term x1 = u.at_time(x, 1);
  ASSERT_NE(x1, x0);
  ASSERT_EQ(x1, u.at_time(x, 1));
}

TEST_P(UnrollerUnitTests, GetTime)
{
  RelationalTransitionSystem rts(s);
  counter_system(rts, rts.make_term(10, bvsort));
  Term x = rts.named_terms().at("x");

  Unroller u(rts);
  Term x0 = u.at_time(x, 0);
  EXPECT_EQ(0, u.get_var_time(x0));

  Term x4 = u.at_time(x, 4);
  EXPECT_EQ(4, u.get_var_time(x4));

  Term x2 = u.at_time(x, 2);
  EXPECT_EQ(2, u.get_var_time(x2));

  // make sure that x4 still works a second time
  Term x4_2 = u.at_time(x, 4);
  EXPECT_EQ(4, u.get_var_time(x4_2));
  EXPECT_EQ(x4, x4_2);

  Term x1 = u.at_time(x, 1);

  Term x0px0 = rts.make_term(BVAdd, x0, x0);
  EXPECT_EQ(0, u.get_curr_time(x0px0));

  Term x1px2 = rts.make_term(BVAdd, x1, x2);
  EXPECT_EQ(1, u.get_curr_time(x1px2));

  Term x1px4 = rts.make_term(BVAdd, x1, x4);
  // cannot get current time for term that could not have been unrolled
  // because difference between times is greater than 1
  // should be <= since only over current state vars, input vars, and next state
  // vars
  EXPECT_THROW(u.get_curr_time(x1px4), PonoException);
}

TEST_P(UnrollerUnitTests, StagedUnrolling)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term y = rts.make_statevar("y", bvsort);

  Unroller u(rts);
  Term xpy = rts.make_term(BVAdd, x, y);
  Term x4py4 = u.at_time(xpy, 4);

  UnorderedTermSet free_vars;
  get_free_symbolic_consts(xpy, free_vars);
  EXPECT_TRUE(free_vars.find(x) != free_vars.end());
  EXPECT_TRUE(free_vars.find(y) != free_vars.end());

  free_vars.clear();
  Term x4 = u.at_time(x, 4);
  Term y4 = u.at_time(y, 4);
  get_free_symbolic_consts(x4py4, free_vars);
  EXPECT_TRUE(free_vars.find(x4) != free_vars.end());
  EXPECT_TRUE(free_vars.find(y4) != free_vars.end());

  // Try staged unrolling
  // meaning unroll some but not all of the variables
  // then unroll the rest all at once

  // x is not unrolled but y is
  Term xpy4 = rts.make_term(BVAdd, x, y4);
  Term x4py4_2 = u.at_time(xpy4, 4);
  EXPECT_EQ(x4py4_2, x4py4);
}

TEST_P(UnrollerUnitTests, Unroller)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);

  Unroller au(rts);
  Term x4 = au.at_time(x, 4);

  // add a new variable after declaring the Unroller
  // this kind of unroller allows this
  Term y = rts.make_statevar("y", bvsort);
  Term xpy = rts.make_term(BVAdd, x, y);

  Term y4 = au.at_time(y, 4);
  // expect it to have actually unrolled
  EXPECT_NE(y4, y);

  Term x4py4 = rts.make_term(BVAdd, x4, y4);
  Term x4py4_2 = au.at_time(xpy, 4);
  EXPECT_EQ(x4py4, x4py4_2);
}

TEST_P(UnrollerUnitTests, FunctionalUnroller)
{
  FunctionalTransitionSystem fts(s);
  counter_system(fts, fts.make_term(10, bvsort));
  Term x = fts.named_terms().at("x");

  // add an input
  Term inp = fts.make_inputvar("inp", bvsort);

  // add another state variable
  Term y = fts.make_statevar("y", bvsort);
  fts.assign_next(y, fts.make_term(BVAdd, y, inp));

  FunctionalUnroller funroller(fts, 0);

  Term x0 = funroller.at_time(x, 0);
  EXPECT_TRUE(x0->is_symbolic_const());
  EXPECT_FALSE(funroller.at_time(x, 1)->is_symbolic_const());

  cout << "Pure functional unrolling to 2 gives: " << funroller.at_time(x, 2)
       << endl;
  cout << "Input unrolled at 2 gives: " << funroller.at_time(inp, 2) << endl;

  EXPECT_TRUE(funroller.at_time(inp, 2)->is_symbolic_const())
      << "inputs always need fresh symbols";

  // only free symbol in pure functional unrolling should be state variables at
  // 0 and inputs since x doesn't have y or inp in the COI, just check that
  UnorderedTermSet free_vars;
  get_free_symbolic_consts(funroller.at_time(x, 3), free_vars);
  EXPECT_EQ(free_vars.size(), 1)
      << "should only be on constant because of pure unrolling";
  EXPECT_TRUE(free_vars.find(x0) != free_vars.end())
      << "expecting one free variable: x@0";

  free_vars.clear();

  // now try with y
  Term y0 = funroller.at_time(y, 0);
  EXPECT_TRUE(y0->is_symbolic_const());
  UnorderedTermSet expected_free_vars({ y0 });
  size_t unroll_length = 4;
  for (size_t i = 0; i < unroll_length; ++i) {
    Term unrolled_inp = funroller.at_time(inp, i);
    EXPECT_TRUE(unrolled_inp->is_symbolic_const());
    expected_free_vars.insert(unrolled_inp);
  }

  Term unrolled_y = funroller.at_time(y, unroll_length);
  get_free_symbolic_consts(unrolled_y, free_vars);
  EXPECT_EQ(free_vars.size(), expected_free_vars.size());
  for (auto fv : expected_free_vars) {
    EXPECT_TRUE(free_vars.find(fv) != free_vars.end()) << "missing " << fv;
  }

  EXPECT_THROW(funroller.at_time(fts.next(x), 4), PonoException)
      << "FunctionalUnroller can't handle next state variables" << endl;

  // check untiming
  // doesn't make tons of sense to untime a functional unrolling
  // nevertheless there are certain times where it's needed
  // just make sure you know what you're doing
  // e.g. in real applications, probably need to substitute
  // in values for input variables BEFORE untiming or you will
  // get a nonsense formula
  Term untimed_y = funroller.untime(unrolled_y);

  free_vars.clear();
  get_free_symbolic_consts(untimed_y, free_vars);
  expected_free_vars.clear();
  expected_free_vars.insert(y);
  expected_free_vars.insert(inp);

  EXPECT_EQ(free_vars.size(), expected_free_vars.size());

  for (auto fv : expected_free_vars) {
    EXPECT_TRUE(free_vars.find(fv) != free_vars.end())
        << "Expected free variable " << fv << " not in untimed term" << endl;
  }
}

TEST_P(UnrollerUnitTests, IntermittentFunctionalUnrolling)
{
  FunctionalTransitionSystem fts(s);
  counter_system(fts, fts.make_term(10, bvsort));
  Term x = fts.named_terms().at("x");

  size_t interval = 4;
  FunctionalUnroller funroller(fts, interval);

  Term x0 = funroller.at_time(x, 0);
  EXPECT_TRUE(x0->is_symbolic_const());

  Term true_ = s->make_term(true);

  for (size_t i = 1; i < interval; ++i) {
    Term unrolled_x = funroller.at_time(x, i);

    UnorderedTermSet free_vars;
    get_free_symbolic_consts(unrolled_x, free_vars);
    EXPECT_EQ(free_vars.size(), 1);
    EXPECT_TRUE(free_vars.find(x0) != free_vars.end());
    EXPECT_EQ(funroller.extra_constraints_at(i), true_);
  }

  Term x_at_interval = funroller.at_time(x, interval);
  EXPECT_TRUE(x_at_interval->is_symbolic_const());

  cout << "Functional unrolling at interval got " << x_at_interval << endl;

  // expected constraint is x@4 = pure functional unrolling at 3 plugged into
  // update function
  Term update = s->substitute(fts.state_updates().at(x),
                              { { x, funroller.at_time(x, 3) } });
  Term expected_eq_constraint = fts.make_term(Equal, x_at_interval, update);

  // not guaranteed to be structurally identical
  // but should be equivalent
  s->assert_formula(s->make_term(Distinct,
                                 expected_eq_constraint,
                                 funroller.extra_constraints_at(interval)));
  Result r = s->check_sat();
  EXPECT_TRUE(r.is_unsat());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedUnrollerUnitTests,
                         UnrollerUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
