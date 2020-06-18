#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

using namespace cosa;
using namespace smt;
using namespace std;

namespace cosa_tests {

class TSUnitTests : public ::testing::Test,
                    public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    boolsort = s->make_sort(BOOL);
    bvsort = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort;
};

TEST(FTSDefaults, DefaultSolverCVC4)
{
  FunctionalTransitionSystem fts;
  ASSERT_TRUE(fts.solver()->get_solver_enum() == CVC4);
}

TEST(RTSDefaults, DefaultSolverCVC4)
{
  RelationalTransitionSystem rts;
  ASSERT_TRUE(rts.solver()->get_solver_enum() == CVC4);
}

TEST_P(TSUnitTests, TransferTransitionSystem)
{
  RelationalTransitionSystem rts;
  SmtSolver solver = rts.solver();

  // Have to use sorts associated with this solver
  // not with the parameterized solver "s"
  Sort bvsort8 = solver->make_sort(BV, 8);
  Sort bsort = solver->make_sort(BOOL);

  Term x = rts.make_statevar("x", bvsort8);
  Term y = rts.make_statevar("y", bvsort8);
  Term incx = rts.make_inputvar("incx", bsort);
  Term incy = rts.make_inputvar("incy", bsort);

  rts.constrain_init(solver->make_term(Equal, x, y));
  rts.assign_next(
      x,
      solver->make_term(
          Ite,
          incx,
          solver->make_term(BVAdd, x, solver->make_term(1, bvsort8)),
          x));

  rts.assign_next(
      y,
      solver->make_term(
          Ite,
          incy,
          solver->make_term(BVAdd, y, solver->make_term(1, bvsort8)),
          y));
  // ((x - y) >= 8) -> !inc_x
  rts.add_constraint(solver->make_term(
      Implies,
      solver->make_term(
          BVUge, solver->make_term(BVSub, x, y), solver->make_term(8, bvsort8)),
      solver->make_term(Not, incx)));
  // ((y - x) >= 8) -> !inc_y
  rts.add_constraint(solver->make_term(
      Implies,
      solver->make_term(
          BVUge, solver->make_term(BVSub, y, x), solver->make_term(8, bvsort8)),
      solver->make_term(Not, incy)));
  RelationalTransitionSystem copied_rts(rts, s);
}

TEST_P(TSUnitTests, FTS_IsFunc)
{
  FunctionalTransitionSystem fts(s);
  ASSERT_TRUE(fts.is_functional());
}

TEST_P(TSUnitTests, RTS_IsFunc)
{
  RelationalTransitionSystem rts(s);
  ASSERT_FALSE(rts.is_functional());
}

TEST_P(TSUnitTests, FTS_Exceptions)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  Term xp1_n = fts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(fts.assign_next(x, xp1_n), CosaException);
}

TEST_P(TSUnitTests, RTS_Exceptions)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term xp1_n = rts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(rts.assign_next(x, xp1_n), CosaException);
  ASSERT_NO_THROW(rts.constrain_trans(s->make_term(Equal, rts.next(x), xp1_n)));
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTSUnitTests,
                         TSUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace cosa_tests
