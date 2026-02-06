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

class MakeTransTotalTests : public ::testing::Test,
                            public ::testing::WithParamInterface<
                                std::tuple<SolverEnum, TSEnum, SortKind>>
{
 protected:
  void SetUp() override
  {
    SolverEnum se = std::get<0>(GetParam());
    TSEnum tse = std::get<1>(GetParam());
    SortKind sk = std::get<2>(GetParam());
    if ((se == BZLA || se == BTOR) && sk == INT) {
      GTEST_SKIP() << "Bitwuzla does not support Integer";
    }
    SmtSolver solver = create_solver(se);
    if (tse == Functional) {
      ts = FunctionalTransitionSystem(solver);
    } else {
      ts = RelationalTransitionSystem(solver);
    }

    sort = (sk == BV) ? ts.make_sort(sk, 8) : ts.make_sort(sk);
    x = ts.make_statevar("x", sort);
    ts.assign_next(
        x, ts.make_term((sk == BV) ? BVAdd : Plus, x, ts.make_term(2, sort)));
    ts.set_init(ts.make_term(Equal, x, ts.make_term(0, sort)));
    ts.add_constraint(
        ts.make_term((sk == BV) ? BVUlt : Lt, x, ts.make_term(10, sort)));
  }
  TransitionSystem ts;
  Sort sort;
  Term x;
};

TEST_P(MakeTransTotalTests, CounterTrue)
{
  Term prop = ts.make_term(
      Equal,
      ts.make_term(
          sort->get_sort_kind() == BV ? BVUrem : Mod, x, ts.make_term(2, sort)),
      ts.make_term(0, sort));
  make_trans_total(ts, prop);
  ASSERT_TRUE(ts.is_right_total());
  ASSERT_TRUE(ts.constraints().empty());
  KInduction kind(SafetyProperty{ ts.solver(), prop }, ts, ts.solver());
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, TRUE);
}

TEST_P(MakeTransTotalTests, CounterFalse)
{
  Term prop = ts.make_term(Not, ts.make_term(Equal, x, ts.make_term(8, sort)));
  make_trans_total(ts, prop);
  ASSERT_TRUE(ts.is_right_total());
  ASSERT_TRUE(ts.constraints().empty());
  KInduction kind(SafetyProperty{ ts.solver(), prop }, ts, ts.solver());
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, FALSE);
  ASSERT_EQ(kind.witness_length(), 4);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedMakeTransTotalTests,
    MakeTransTotalTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(std::vector<TSEnum>{ Functional,
                                                            Relational }),
                     testing::ValuesIn(std::vector<SortKind>{ BV, INT })));

}  // namespace pono_tests
