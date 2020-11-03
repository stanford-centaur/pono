#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/term_analysis.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class TermAnalysisUnitTests : public ::testing::Test,
                              public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    boolsort = s->make_sort(BOOL);
    bvsort = s->make_sort(BV, 8);
    funsort = s->make_sort(FUNCTION, { bvsort, bvsort, bvsort });
    relsort = s->make_sort(FUNCTION, { bvsort, bvsort, boolsort });
  }
  SmtSolver s;
  Sort boolsort, bvsort, funsort, relsort;
};

TEST_P(TermAnalysisUnitTests, GetPredicates)
{
  if (s->get_solver_enum() == BTOR) {
    std::cout << "SKIPPING get predicates test for BTOR" << std::endl;
    std::cout << "smt-switch btor backend cannot get sort of a UF" << std::endl;
    std::cout << "which is required for get_predicates" << std::endl;
    return;
  }

  Term x = s->make_symbol("x", bvsort);
  Term y = s->make_symbol("y", bvsort);
  Term b = s->make_symbol("b", boolsort);
  Term f = s->make_symbol("f", funsort);
  Term R = s->make_symbol("R", relsort);

  Term xpy = s->make_term(BVAdd, x, y);
  Term fxy = s->make_term(Apply, f, x, y);
  Term Rxy = s->make_term(Apply, R, x, y);
  Term bRxy = s->make_term(And, b, Rxy);
  Term xpy_eq_fxy = s->make_term(Equal, xpy, fxy);
  Term formula = s->make_term(Implies, bRxy, xpy_eq_fxy);

  UnorderedTermSet expected_preds({ Rxy, xpy_eq_fxy });
  UnorderedTermSet preds;
  // by default won't include boolean symbols
  get_predicates(formula, boolsort, preds);
  EXPECT_EQ(preds.size(), expected_preds.size());
  // check they're the same
  for (auto p : preds) {
    EXPECT_TRUE(expected_preds.find(p) != expected_preds.end());
  }

  // boolean symbol b should not have been included
  EXPECT_TRUE(preds.find(b) == preds.end());

  // now try adding the boolean symbol
  get_predicates(formula, boolsort, preds, true);
  // added one predicate -- b
  EXPECT_TRUE(preds.size() == expected_preds.size() + 1);
  EXPECT_TRUE(preds.find(b) != preds.end());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedTermAnalysisUnitTests,
                         TermAnalysisUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
