#include "core/fts.h"
#include "gtest/gtest.h"
#include "modifiers/static_coi.h"
#include "smt/available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class CoiUnitTests : public ::testing::Test,
                     public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort8;
};

TEST_P(CoiUnitTests, SimpleCoiTest)
{
  FunctionalTransitionSystem fts(s);

  Term a = fts.make_inputvar("a", bvsort8);
  Term b = fts.make_inputvar("b", bvsort8);
  Term c = fts.make_inputvar("c", bvsort8);
  Term d = fts.make_inputvar("d", boolsort);

  Term regres = fts.make_statevar("regres", bvsort8);
  Term counter = fts.make_statevar("counter", bvsort8);

  Term zero = fts.make_term(0, bvsort8);
  Term one = fts.make_term(1, bvsort8);
  fts.constrain_init(fts.make_term(Equal, counter, zero));
  fts.assign_next(counter, fts.make_term(BVAdd, counter, one));
  fts.assign_next(regres, fts.make_term(BVAdd, a, b));
  fts.add_constraint(fts.make_term(Equal, d, fts.make_term(BVUle, a, b)));

  Term out = fts.make_term(BVSub, regres, one);
  fts.name_term("out", out);

  StaticConeOfInfluence coi(fts, { out });

  const UnorderedTermSet & statevars = fts.statevars();
  const UnorderedTermSet & inputvars = fts.inputvars();
  const unordered_map<string, Term> & named_terms = fts.named_terms();
  EXPECT_EQ(inputvars.size(), 3);
  EXPECT_TRUE(inputvars.find(a) != inputvars.end());
  EXPECT_TRUE(inputvars.find(b) != inputvars.end());
  EXPECT_TRUE(inputvars.find(c) == inputvars.end());
  // d is included because it is part of a constraint
  EXPECT_TRUE(inputvars.find(d) != inputvars.end());
  EXPECT_TRUE(named_terms.find("a") != named_terms.end());
  EXPECT_TRUE(named_terms.find("c") == named_terms.end());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedCoiUnitTests,
                         CoiUnitTests,
                         testing::ValuesIn(available_solver_enums()));
}  // namespace pono_tests
