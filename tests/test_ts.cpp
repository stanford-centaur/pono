#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/prop.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class TSUnitTests : public ::testing::Test,
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

TEST_P(TSUnitTests, FTS_IsFunc)
{
  FunctionalTransitionSystem fts(s);
  ASSERT_TRUE(fts.is_functional());

  // state variables without state updates
  // will make the system non-deterministic
  Term x = fts.make_statevar("x", bvsort);
  ASSERT_FALSE(fts.is_deterministic());
  ASSERT_TRUE(fts.is_functional());

  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_TRUE(fts.is_functional());
  ASSERT_TRUE(fts.is_deterministic());

  fts.add_constraint(fts.make_term(BVUge, x, s->make_term(2, bvsort)));
  // any kind of constrains makes the system non-deterministic
  // TODO need to improve names here
  ASSERT_FALSE(fts.is_deterministic());

  Term y = fts.make_statevar("y", bvsort);
  fts.assign_next(y, y);
  ASSERT_TRUE(fts.is_functional());
  // still can't be deterministic because of the constraint
  ASSERT_FALSE(fts.is_deterministic());

  TransitionSystem ts_copy = fts;
  ASSERT_TRUE(ts_copy.is_functional());
}

TEST_P(TSUnitTests, RTS_IsFunc)
{
  RelationalTransitionSystem rts(s);
  ASSERT_FALSE(rts.is_functional());

  // state variables without state updates
  // will make the system non-functional
  Term x = rts.make_statevar("x", bvsort);
  ASSERT_FALSE(rts.is_functional());

  rts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  // Relational transition system is still not functional
  ASSERT_FALSE(rts.is_functional());
  // cannot guarantee determinism if relational
  ASSERT_FALSE(rts.is_deterministic());

  TransitionSystem ts_copy = rts;
  ASSERT_FALSE(ts_copy.is_functional());
}

TEST_P(TSUnitTests, FTS_Exceptions)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  Term xp1_n = fts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(fts.assign_next(x, xp1_n), PonoException);
}

TEST_P(TSUnitTests, RTS_Exceptions)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term xp1_n = rts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(rts.assign_next(x, xp1_n), PonoException);
  ASSERT_NO_THROW(rts.constrain_trans(s->make_term(Equal, rts.next(x), xp1_n)));
}

TEST_P(TSUnitTests, RTS_Copy)
{
  RelationalTransitionSystem rts(s);

  RelationalTransitionSystem rts2 = rts;
  TransitionSystem ts = rts;
}

TEST_P(TSUnitTests, Prop_Copy)
{
  RelationalTransitionSystem rts(s);
  Property p(rts, s->make_term(true));

  Property p2 = p;
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTSUnitTests,
                         TSUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
