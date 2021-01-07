#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/term_analysis.h"
#include "utils/term_walkers.h"

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

TEST_P(TermAnalysisUnitTests, GetPredicatesBasic)
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
  get_predicates(s, formula, preds);
  EXPECT_EQ(preds.size(), expected_preds.size());
  // check they're the same
  for (auto p : preds) {
    EXPECT_TRUE(expected_preds.find(p) != expected_preds.end());
  }

  // boolean symbol b should not have been included
  EXPECT_TRUE(preds.find(b) == preds.end());

  // now try adding the boolean symbol
  get_predicates(s, formula, preds, true);
  // added one predicate -- b
  EXPECT_TRUE(preds.size() == expected_preds.size() + 1);
  EXPECT_TRUE(preds.find(b) != preds.end());
}

TEST_P(TermAnalysisUnitTests, GetPredicatesIte)
{
  Term x = s->make_symbol("x", bvsort);
  Term y = s->make_symbol("y", bvsort);
  Term z = s->make_symbol("z", bvsort);
  Term nextval = s->make_symbol("nextval", bvsort);

  Term xley = s->make_term(BVUle, x, y);
  Term ylt8 = s->make_term(BVUlt, y, s->make_term(8, bvsort));
  Term yltz = s->make_term(BVUlt, y, z);
  Term update = s->make_term(
      Ite, s->make_term(And, xley, ylt8), x, s->make_term(Ite, yltz, y, z));
  Term formula = s->make_term(Equal, nextval, update);

  // expecting the normal predicates + removing the ITEs
  // thus for each ITE where a = ite(p, b, c) we get
  // p as a predicate plus a = b and a = c
  // and we don't get any ITEs in our predicates
  UnorderedTermSet eq_preds({ s->make_term(Equal, nextval, x),
                              s->make_term(Equal, nextval, y),
                              s->make_term(Equal, nextval, z) });
  UnorderedTermSet expected_preds({ xley, ylt8, yltz });
  for (auto ep : eq_preds) {
    expected_preds.insert(ep);
  }
  UnorderedTermSet preds;
  // use option to split ITEs
  bool include_symbols = false;
  bool search_subterms = false;
  bool split_ites = true;
  get_predicates(
      s, formula, preds, include_symbols, search_subterms, split_ites);

  EXPECT_EQ(preds.size(), expected_preds.size());

  // due to rewriting, they may not be exactly the same predicates
  // e.g. x <= y gets rewritten to not(y < x) and then the predicate is y < x
  // However, the equality predicates should not be rewritten so we can check
  // those
  for (auto ep : eq_preds) {
    // expected equal predicate should be in the found predicates
    EXPECT_TRUE(preds.find(ep) != preds.end());
  }

  // Now, also check that there are no ITEs in the found predicates
  TermOpCollector opfinder(s);
  for (auto p : preds) {
    UnorderedTermSet out;
    opfinder.find_matching_terms(p, { Ite }, out);
    EXPECT_TRUE(out.empty());  // expecting no ITEs in the found predicates
  }
}

INSTANTIATE_TEST_SUITE_P(ParameterizedTermAnalysisUnitTests,
                         TermAnalysisUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
