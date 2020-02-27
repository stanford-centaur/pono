#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/unroller.h"
#include "utils/exceptions.h"

#include "smt-switch/boolector_factory.h"
#include "smt-switch/smt.h"

#include "available_solvers.h"

using namespace cosa;
using namespace smt;
using namespace std;

namespace cosa_tests {

class UnitTests : public ::testing::Test,
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

TEST_P(UnitTests, FTS_Unroll)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_state("x", bvsort);
  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));

  Unroller u(fts, s);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
}

TEST_P(UnitTests, FTS_Exceptions)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_state("x", bvsort);
  Term xp1_n = fts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  ASSERT_THROW(fts.assign_next(x, xp1_n), CosaException);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverUnitTests,
                         UnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace cosa_tests
