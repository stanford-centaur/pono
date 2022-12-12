#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/rts.h"
#include "engines/ic3.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/ts_analysis.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class IC3UnitTests : public ::testing::Test,
                     public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver_for(GetParam(), IC3_BOOL, false);
    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort8;
};

TEST_P(IC3UnitTests, SimpleSystemSafe)
{
  RelationalTransitionSystem rts(s);
  Term s1 = rts.make_statevar("s1", boolsort);
  Term s2 = rts.make_statevar("s2", boolsort);

  // INIT !s1 & !s2
  rts.constrain_init(s->make_term(Not, s1));
  rts.constrain_init(s->make_term(Not, s2));

  // TRANS next(s1) = (s1 | s2)
  // TRANS next(s2) = s2
  rts.assign_next(s1, s->make_term(Or, s1, s2));
  rts.assign_next(s2, s2);

  Property p(s, s->make_term(Not, s1));

  IC3 ic3(p, rts, s);
  ProverResult r = ic3.prove();
  ASSERT_EQ(r, TRUE);

  // get the invariant
  Term invar = ic3.invar();
  ASSERT_TRUE(check_invar(rts, p.prop(), invar));
}

TEST_P(IC3UnitTests, SimpleSystemUnsafe)
{
  FunctionalTransitionSystem fts(s);
  Term s1 = fts.make_statevar("s1", boolsort);
  Term s2 = fts.make_statevar("s2", boolsort);

  // INIT !s1 & s2
  fts.constrain_init(s->make_term(Not, s1));
  fts.constrain_init(s2);

  // TRANS next(s1) = (s1 | s2)
  // TRANS next(s2) = s2
  fts.assign_next(s1, s->make_term(Or, s1, s2));
  fts.assign_next(s2, s2);

  Property p(s, s->make_term(Not, s1));

  IC3 ic3(p, fts, s);
  ProverResult r = ic3.prove();
  ASSERT_EQ(r, FALSE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverIC3UnitTests,
    IC3UnitTests,
    testing::ValuesIn(available_solver_enums()));
}  // namespace pono_tests
