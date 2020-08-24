#include <utility>
#include <vector>

#include "gtest/gtest.h"

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

INSTANTIATE_TEST_SUITE_P(ParameterizedUnrollerUnitTests,
                         UnrollerUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
