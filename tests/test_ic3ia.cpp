#include <utility>
#include <vector>

#include "available_solvers.h"
#include "core/fts.h"
#include "engines/ic3ia.h"
#include "gtest/gtest.h"
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
    s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    s->set_opt("produce-unsat-cores", "true");
    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort8;
};

TEST_P(IC3UnitTests, SimpleSystemSafe)
{
  FunctionalTransitionSystem fts(s);
  Term s1 = fts.make_statevar("s1", boolsort);
  Term s2 = fts.make_statevar("s2", boolsort);

  // INIT !s1 & !s2
  fts.constrain_init(s->make_term(Not, s1));
  fts.constrain_init(s->make_term(Not, s2));

  // TRANS next(s1) = (s1 | s2)
  // TRANS next(s2) = s2
  fts.assign_next(s1, s->make_term(Or, s1, s2));
  fts.assign_next(s2, s2);

  Property p(fts, s->make_term(Not, s1));

  IC3IA ic3ia(p, s);
  ProverResult r = ic3ia.prove();
  ASSERT_EQ(r, TRUE);

  // get the invariant
  Term invar = ic3ia.invar();
  ASSERT_TRUE(check_invar(fts, p.prop(), invar));
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

  Property p(fts, s->make_term(Not, s1));

  IC3IA ic3ia(p, s);
  ProverResult r = ic3ia.prove();
  ASSERT_EQ(r, FALSE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverIC3UnitTests,
    IC3UnitTests,
    // only using MathSAT for now, but could be more general in the future
    testing::ValuesIn({ MSAT }));
}  // namespace pono_tests
