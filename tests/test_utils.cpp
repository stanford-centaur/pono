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
#include "utils/term_analysis.h"
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
  Property prop(fts.solver(), prop_term);

  SolverEnum se = get<0>(GetParam());
  Engine eng = get<1>(GetParam());

  if (eng == INTERP && se != MSAT) {
    // skip interpolation unless the solver is MathSAT
    return;
  }

  SmtSolver s = create_solver(se);
  s->set_opt("produce-unsat-assumptions", "true");
  std::shared_ptr<Prover> prover = make_prover(eng, prop, fts, s);
  ProverResult r = prover->check_until(9);

  ASSERT_EQ(r, FALSE);
}

TEST_P(UtilsUnitTests, RemoveItes)
{
  Term a = s->make_symbol("a", boolsort);
  Term b = s->make_symbol("b", boolsort);
  Term x = s->make_symbol("x", bvsort);
  Term y = s->make_symbol("y", bvsort);
  Term one = s->make_term(1, bvsort);

  Term xp1 = s->make_term(BVAdd, x, one);
  Term yp1 = s->make_term(BVAdd, y, one);
  Term xpy = s->make_term(BVAdd, x, y);

  TermVec conditions(
      { s->make_term(BVUlt, xp1, s->make_term(10, bvsort)),
        s->make_term(BVUge, yp1, s->make_term(40, bvsort)),
        s->make_term(
            And,
            b,
            s->make_term(
                And, a, s->make_term(BVUge, xpy, s->make_term(20, bvsort)))) });

  ASSERT_EQ(conditions.size(), 3);

  Term ite_1 = s->make_term(Ite, conditions[0], xp1, x);
  Term ite_2 = s->make_term(Ite, conditions[1], y, yp1);
  Term ite = s->make_term(Ite, conditions[2], ite_1, ite_2);

  Term top = s->make_term(BVSub, s->make_term(BVSub, ite, one), one);

  s->assert_formula(s->make_term(Or, a, s->make_term(Not, b)));

  Result r = s->check_sat();
  EXPECT_TRUE(r.is_sat());

  TermVec replaced = remove_ites_under_model(s, { top });
  ASSERT_EQ(replaced.size(), 1);
  Term symbolic_return_val = replaced[0];
  EXPECT_FALSE(symbolic_return_val->is_value());
  EXPECT_EQ(s->get_value(top), s->get_value(symbolic_return_val));

  Term s_true = s->make_term(true);
  TermVec assertions;
  assertions.reserve(conditions.size());
  for (const auto & c : conditions) {
    ASSERT_EQ(c->get_sort(), boolsort);
    if (s->get_value(c) == s_true) {
      assertions.push_back(c);
    } else {
      assertions.push_back(s->make_term(Not, c));
    }
  }

  for (const auto & a : assertions) {
    s->assert_formula(a);
  }

  // check that these terms are equivalent under this path
  s->assert_formula(s->make_term(Distinct, top, symbolic_return_val));
  r = s->check_sat();
  EXPECT_TRUE(r.is_unsat());
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
