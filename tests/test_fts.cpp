#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/unroller.h"
#include "gtest/gtest.h"
#include "smt-switch/boolector_factory.h"
#include "smt-switch/smt.h"

using namespace cosa;
using namespace smt;
using namespace std;

namespace cosa_tests {

class UnitTests : public ::testing::Test
{
 protected:
  void SetUp() override
  {
    s = BoolectorSolverFactory::create();
    bvsort = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort bvsort;
};

TEST_F(UnitTests, FTS)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_state("x", bvsort);
  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));

  Unroller u(fts, s);
  Term x0 = u.at_time(x, 0);
  ASSERT_EQ(x0, u.at_time(x, 0));
}

}  // namespace cosa_tests
