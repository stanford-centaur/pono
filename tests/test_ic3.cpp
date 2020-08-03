#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "engines/mbic3.h"

#include "available_solvers.h"

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
    s = available_solvers().at(GetParam())(false);
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
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

  ModelBasedIC3 mbic3(p, s);
  ProverResult r = mbic3.prove();
  ASSERT_EQ(r, TRUE);
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

  ModelBasedIC3 mbic3(p, s);
  ProverResult r = mbic3.prove();
  ASSERT_EQ(r, FALSE);
}

// TODO: Test over all solvers once boolector's unsat core is fixed

vector<SolverEnum> get_solver_enums()
{
  vector<SolverEnum> solver_enums;
  for (auto se : available_solver_enums()) {
    if (se != BTOR) {
      solver_enums.push_back(se);
    }
  }
  return solver_enums;
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverIC3UnitTests,
                         IC3UnitTests,
                         testing::ValuesIn(get_solver_enums()));
}  // namespace pono_tests
