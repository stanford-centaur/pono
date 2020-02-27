#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/unroller.h"
#include "gtest/gtest.h"
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

TEST_P(UnitTests, FTS)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_state("x", bvsort);
  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));

  Unroller u(fts, s);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverUnitTests,
                         UnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace cosa_tests
