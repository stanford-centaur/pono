#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/exceptions.h"
#include "utils/make_provers.h"
#include "utils/term_walkers.h"
#include "utils/ts_analysis.h"

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

class UtilsEngineUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, Engine>>
{
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

TEST_P(UtilsUnitTests, CheckInvarTrue)
{
  RelationalTransitionSystem rts(s);
  counter_system(rts, rts.make_term(10, bvsort));
  Term x = rts.named_terms().at("x");

  Term prop = rts.make_term(BVUle, x, rts.make_term(10, bvsort));
  // property is inductive
  EXPECT_TRUE(check_invar(rts, prop, prop));
}

TEST_P(UtilsUnitTests, CheckInvarFalse)
{
  RelationalTransitionSystem rts(s);
  counter_system(rts, rts.make_term(11, bvsort));
  Term x = rts.named_terms().at("x");

  Term prop = rts.make_term(BVUle, x, rts.make_term(10, bvsort));
  EXPECT_FALSE(check_invar(rts, prop, prop));
}

TEST_P(UtilsUnitTests, CheckInvarFalseNonState)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  rts.constrain_init(rts.make_term(Equal, x, rts.make_term(0, bvsort)));
  Term inp = rts.make_inputvar("inp", bvsort);
  // x' = x + inp < 10 ? x + inp : 0
  Term xpinp = rts.make_term(BVAdd, x, inp);
  rts.assign_next(
      x,
      rts.make_term(Ite,
                    rts.make_term(BVUlt, xpinp, rts.make_term(10, bvsort)),
                    xpinp,
                    rts.make_term(0, bvsort)));

  Term prop = rts.make_term(BVUle, x, rts.make_term(10, bvsort));

  // invariant includes inputs -- which is not allowed in an invariant
  Term invar = rts.make_term(BVUle, xpinp, rts.make_term(10, bvsort));

  EXPECT_FALSE(check_invar(rts, prop, invar));
}

TEST_P(UtilsEngineUnitTests, MakeProver)
{
  // use default solver
  FunctionalTransitionSystem fts;
  Sort bvsort8 = fts.make_sort(BV, 8);
  Sort boolsort = fts.make_sort(BOOL);
  Term one = fts.make_term(1, bvsort8);
  Term eight = fts.make_term(8, bvsort8);
  Term x = fts.make_statevar("x", bvsort8);

  fts.set_init(fts.make_term(Equal, x, fts.make_term(0, bvsort8)));
  fts.assign_next(x, fts.make_term(BVAdd, x, one));

  Term prop_term = fts.make_term(BVUlt, x, eight);
  Property prop(fts, prop_term);

  SolverEnum se = get<0>(GetParam());
  Engine eng = get<1>(GetParam());

  if (eng == INTERP && se != MSAT) {
    // skip interpolation unless the solver is MathSAT
    return;
  }

  SmtSolver s = create_solver(se);
  s->set_opt("incremental", "true");
  s->set_opt("produce-models", "true");

  ProverResult r;
  if (eng == INTERP && se == MSAT) {
    SmtSolver interp_s = create_interpolating_solver(MSAT_INTERPOLATOR);
    std::shared_ptr<Prover> prover = make_prover(eng, prop, s, interp_s);
    r = prover->check_until(9);
  } else {
    std::shared_ptr<Prover> prover = make_prover(eng, prop, s);
    r = prover->check_until(9);
  }
  ASSERT_EQ(r, FALSE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedUtilsUnitTests,
                         UtilsUnitTests,
                         testing::ValuesIn(available_solver_enums()));

INSTANTIATE_TEST_SUITE_P(
    ParameterizedUtilsEngineUnitTests,
    UtilsEngineUnitTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(all_engines())));

}  // namespace pono_tests
