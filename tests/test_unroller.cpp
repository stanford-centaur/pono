#include <utility>
#include <vector>

#include "smt-switch/utils.h"

#include "gtest/gtest.h"

#include "core/adaptive_unroller.h"
#include "core/fts.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

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
  Term x = fts.make_statevar("x", bvsort);
  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));

  Unroller u(fts, s);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
}

TEST_P(UnrollerUnitTests, RTS_Unroll)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  rts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));

  Unroller u(rts, s);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
  Term x1 = u.at_time(x, 1);
  ASSERT_NE(x1, x0);
  ASSERT_EQ(x1, u.at_time(x, 1));
}

TEST_P(UnrollerUnitTests, GetTime)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  rts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));

  Unroller u(rts, s);
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

  Unroller u(rts, s);
  Term xpy = rts.make_term(BVAdd, x, y);
  Term x4py4 = u.at_time(xpy, 4);

  TermVec free_vars_vec;
  UnorderedTermSet free_vars;
  get_free_symbolic_consts(xpy, free_vars_vec);
  free_vars.insert(free_vars_vec.begin(), free_vars_vec.end());
  EXPECT_TRUE(free_vars.find(x) != free_vars.end());
  EXPECT_TRUE(free_vars.find(y) != free_vars.end());

  free_vars_vec.clear();
  free_vars.clear();
  Term x4 = u.at_time(x, 4);
  Term y4 = u.at_time(y, 4);
  get_free_symbolic_consts(x4py4, free_vars_vec);
  free_vars.insert(free_vars_vec.begin(), free_vars_vec.end());
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

TEST_P(UnrollerUnitTests, AdaptiveUnroller)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);

  AdaptiveUnroller au(rts, s);
  Term x4 = au.at_time(x, 4);

  // add a new variable after declaring the AdaptiveUnroller
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

INSTANTIATE_TEST_SUITE_P(ParameterizedUnrollerUnitTests,
                         UnrollerUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
