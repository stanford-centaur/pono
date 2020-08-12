#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/rts.h"
#include "core/unroller.h"
#include "engines/kinduction.h"
#include "utils/exceptions.h"
#include "utils/term_walkers.h"

#include "available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class UtilsUnitTests : public ::testing::Test,
                       public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    boolsort = s->make_sort(BOOL);
    bvsort = s->make_sort(BV, 8);
    funsort = s->make_sort(FUNCTION, { bvsort, boolsort });
  }
  SmtSolver s;
  Sort boolsort, bvsort, funsort;
};

TEST_P(UtilsUnitTests, FindApply)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term y = rts.make_statevar("y", bvsort);
  Term xm1 = rts.make_term(BVSub, x, rts.make_term(1, bvsort));
  Term f = rts.solver()->make_symbol("f", funsort);
  Term g = rts.solver()->make_symbol("g", funsort);

  Term fx = rts.make_term(Apply, f, x);
  Term gxm1 = rts.make_term(Apply, g, xm1);
  Term fy = rts.make_term(Apply, f, y);
  Term xpy = rts.make_term(BVAdd, x, y);
  Term gxpy = rts.make_term(Apply, g, xpy);

  Term term = rts.make_term(Implies, fx, gxm1);
  term = rts.make_term(Or, term, gxpy);
  term = rts.make_term(And, term, fy);

  TermOpCollector toc(s);
  // find terms with op Apply and BVSub
  UnorderedTermSet matching_terms;
  toc.find_matching_terms(term, { Apply }, matching_terms);
  EXPECT_EQ(matching_terms.size(), 4);
  EXPECT_TRUE(matching_terms.find(fx) != matching_terms.end());
  EXPECT_TRUE(matching_terms.find(fy) != matching_terms.end());
  EXPECT_TRUE(matching_terms.find(gxm1) != matching_terms.end());
  EXPECT_TRUE(matching_terms.find(gxpy) != matching_terms.end());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedUtilsUnitTests,
                         UtilsUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
