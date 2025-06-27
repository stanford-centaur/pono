#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/term_walkers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class WalkersUnitTests : public ::testing::Test,
                         public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    boolsort = s->make_sort(BOOL);
    bvsort4 = s->make_sort(BV, 8);
    bvsort8 = s->make_sort(BV, 8);
    funsort = s->make_sort(FUNCTION, { bvsort8, bvsort8, bvsort4 });
    relsort = s->make_sort(FUNCTION, { bvsort4, bvsort8, boolsort });

    x = s->make_symbol("x", bvsort8);
    y = s->make_symbol("y", bvsort8);
    z = s->make_symbol("z", bvsort4);
    f = s->make_symbol("f", funsort);
    r = s->make_symbol("r", relsort);
  }
  SmtSolver s;
  Sort boolsort, bvsort4, bvsort8, funsort, relsort;
  Term x, y, z, f, r;
};

TEST_P(WalkersUnitTests, SubTermCollectorBasic)
{
  Term f_xy = s->make_term(Apply, f, x, y);
  Term r_fxy_y = s->make_term(Apply, r, f_xy, y);
  Term ite_term = s->make_term(Ite, r_fxy_y, f_xy, z);
  Term sum = s->make_term(BVAdd, ite_term, z);

  // manually collect all subterms
  // don't worry about DAG traversal because it's small
  UnorderedTermSet all_subterms;
  TermVec to_visit({ sum });
  while (to_visit.size()) {
    Term t = to_visit.back();
    to_visit.pop_back();
    all_subterms.insert(t);
    for (auto tt : t) {
      to_visit.push_back(tt);
    }
  }

  // default behavior is to exclude functions and ITEs
  SubTermCollector stc(s);
  stc.collect_subterms(sum);

  const unordered_map<Sort, UnorderedTermSet> & stc_subterms =
      stc.get_subterms();
  EXPECT_TRUE(stc_subterms.find(boolsort) == stc_subterms.end());
  EXPECT_TRUE(stc_subterms.find(funsort) == stc_subterms.end());
  EXPECT_TRUE(stc_subterms.find(relsort) == stc_subterms.end());

  Sort sort;
  for (auto term : all_subterms) {
    sort = term->get_sort();
    if (term->get_op() == Ite) {
      EXPECT_TRUE(stc_subterms.find(sort) != stc_subterms.end())
          << "expecting to have terms of the same sort because of ITE children"
          << endl;
      EXPECT_TRUE(stc_subterms.at(sort).find(term)
                  == stc_subterms.at(sort).end())
          << "expecting no ITE terms to be saved" << endl;
    } else if (sort->get_sort_kind() != FUNCTION && sort != boolsort) {
      EXPECT_TRUE(stc_subterms.find(sort) != stc_subterms.end())
          << "Missing sort " << sort << endl;
      EXPECT_TRUE(stc_subterms.at(sort).find(term)
                  != stc_subterms.at(sort).end());
    }
  }

  // now try collecting with boolean terms, functions and ITEs also
  SubTermCollector stcfun(s, false, false, false);
  stcfun.collect_subterms(sum);

  const unordered_map<Sort, UnorderedTermSet> & stcfun_subterms =
      stcfun.get_subterms();
  EXPECT_TRUE(stcfun_subterms.find(funsort) != stcfun_subterms.end());
  EXPECT_TRUE(stcfun_subterms.find(relsort) != stcfun_subterms.end());

  for (auto term : all_subterms) {
    sort = term->get_sort();
    if (sort != boolsort) {
      EXPECT_TRUE(stcfun_subterms.find(sort) != stcfun_subterms.end())
          << "Missing sort " << sort << endl;
      EXPECT_TRUE(stcfun_subterms.at(sort).find(term)
                  != stcfun_subterms.at(sort).end());
    } else {
      EXPECT_TRUE(stcfun.get_predicates().find(term)
                  != stcfun.get_predicates().end())
          << "Missing predicate " << term << endl;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(ParameterizedWalkersUnitTests,
                         WalkersUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
