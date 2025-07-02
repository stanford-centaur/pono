#include "core/fts.h"
#include "core/prop.h"
#include "core/rts.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"

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
  EXPECT_TRUE(fts.is_functional());

  // state variables without state updates
  // will make the system non-deterministic
  Term x = fts.make_statevar("x", bvsort);
  EXPECT_FALSE(fts.is_deterministic());
  EXPECT_TRUE(fts.is_functional());
  EXPECT_EQ(fts.statevars_with_no_update(), UnorderedTermSet({ x }));

  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  EXPECT_TRUE(fts.is_functional());
  EXPECT_TRUE(fts.is_deterministic());
  EXPECT_EQ(fts.statevars_with_no_update().size(), 0);

  fts.add_constraint(fts.make_term(BVUge, x, s->make_term(2, bvsort)));
  // any kind of constrains makes the system non-deterministic
  // TODO need to improve names here
  EXPECT_FALSE(fts.is_deterministic());

  Term y = fts.make_statevar("y", bvsort);
  fts.assign_next(y, y);
  EXPECT_TRUE(fts.is_functional());
  // still can't be deterministic because of the constraint
  EXPECT_FALSE(fts.is_deterministic());
  EXPECT_EQ(fts.statevars_with_no_update().size(), 0);

  Term ynext = fts.make_statevar("y.next", bvsort);
  EXPECT_EQ(fts.statevars_with_no_update(), UnorderedTermSet({ ynext }));

  TransitionSystem ts_copy = fts;
  EXPECT_EQ(fts.is_functional(), ts_copy.is_functional());
  EXPECT_EQ(fts.is_deterministic(), ts_copy.is_deterministic());
  EXPECT_EQ(fts.statevars_with_no_update(), ts_copy.statevars_with_no_update());
  EXPECT_EQ(ts_copy, fts);
  ts_copy.set_init(ts_copy.make_term(Equal, x, ts_copy.make_term(1, bvsort)));
  EXPECT_NE(ts_copy, fts);
}

TEST_P(TSUnitTests, RTS_IsFunc)
{
  RelationalTransitionSystem rts(s);
  EXPECT_FALSE(rts.is_functional());

  // state variables without state updates
  // will make the system non-functional
  Term x = rts.make_statevar("x", bvsort);
  EXPECT_FALSE(rts.is_functional());

  rts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  // Relational transition system is still not functional
  EXPECT_FALSE(rts.is_functional());
  // cannot guarantee determinism if relational
  EXPECT_FALSE(rts.is_deterministic());

  TransitionSystem ts_copy = rts;
  EXPECT_EQ(rts.is_functional(), ts_copy.is_functional());
  EXPECT_EQ(rts.is_deterministic(), ts_copy.is_deterministic());
  EXPECT_EQ(ts_copy, rts);
  ts_copy.set_init(ts_copy.make_term(Equal, x, ts_copy.make_term(1, bvsort)));
  EXPECT_NE(ts_copy, rts);
}

TEST_P(TSUnitTests, FTS_Exceptions)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  Term xp1_n = fts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  EXPECT_THROW(fts.assign_next(x, xp1_n), PonoException);
}

TEST_P(TSUnitTests, RTS_Exceptions)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term xp1_n = rts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  EXPECT_THROW(rts.assign_next(x, xp1_n), PonoException);
  EXPECT_NO_THROW(rts.constrain_trans(s->make_term(Equal, rts.next(x), xp1_n)));
}

TEST_P(TSUnitTests, FTS_DefaultCopy)
{
  FunctionalTransitionSystem fts;
  EXPECT_NO_THROW(fts = FunctionalTransitionSystem(s));
  // make sure the terms are not null
  EXPECT_TRUE(fts.init());
  EXPECT_TRUE(fts.trans());
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
  SafetyProperty p(s, s->make_term(true));

  SafetyProperty p2 = p;
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTSUnitTests,
                         TSUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
