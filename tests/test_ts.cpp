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
    s = available_solvers().at(GetParam())();
    bvsort = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort bvsort;
};

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
  Term x = fts.make_state("x", bvsort);
  Term xp1_n = fts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(fts.assign_next(x, xp1_n), CosaException);
}

TEST_P(TSUnitTests, RTS_Exceptions)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_state("x", bvsort);
  Term xp1_n = rts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(rts.assign_next(x, xp1_n), CosaException);
  ASSERT_NO_THROW(rts.constrain_trans(s->make_term(Equal, rts.next(x), xp1_n)));
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTSUnitTests,
                         TSUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace cosa_tests
