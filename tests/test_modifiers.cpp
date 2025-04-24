#include <utility>

#include "core/fts.h"
#include "core/rts.h"
#include "gtest/gtest.h"
#include "modifiers/history_modifier.h"
#include "modifiers/implicit_predicate_abstractor.h"
#include "modifiers/liveness_to_safety_translator.h"
#include "modifiers/prophecy_modifier.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class ModifierUnitTests : public ::testing::Test,
                          public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    boolsort = s->make_sort(BOOL);
    bvsort = s->make_sort(BV, 8);
    arrsort = s->make_sort(ARRAY, bvsort, bvsort);
  }
  SmtSolver s;
  Sort boolsort, bvsort, arrsort;
};

TEST_P(ModifierUnitTests, HistoryModifier)
{
  FunctionalTransitionSystem fts(s);
  Term max_val = fts.make_term(10, bvsort);
  counter_system(fts, max_val);
  Term x = fts.named_terms().at("x");

  HistoryModifier hm(fts);

  size_t num_state_vars_orig = fts.statevars().size();

  Term trans_1 = fts.trans();
  Term hist_x_10 = hm.get_hist(x, 10);
  Term trans_2 = fts.trans();

  // should have added history variables to the transition relation
  EXPECT_NE(trans_1, trans_2);
  EXPECT_EQ(num_state_vars_orig + 10, fts.statevars().size());

  // shouldn't need to modify system
  // this one already needed to be created
  Term hist_x_2 = hm.get_hist(x, 2);
  Term trans_3 = fts.trans();

  EXPECT_EQ(trans_2, trans_3);
  EXPECT_EQ(num_state_vars_orig + 10, fts.statevars().size());
}

TEST_P(ModifierUnitTests, ProphecyModifierSimple)
{
  FunctionalTransitionSystem fts(s);
  counter_system(fts, fts.make_term(9, bvsort));
  Term x = fts.named_terms().at("x");

  Term prop = fts.make_term(BVUlt, x, fts.make_term(10, bvsort));
  size_t num_statevars_orig = fts.statevars().size();
  EXPECT_EQ(num_statevars_orig, 1);

  ProphecyModifier pm(fts);
  std::pair<Term, Term> p = pm.get_proph(x, 2);
  Term proph_var = p.first;
  Term new_target = p.second;
  // update the property with the prophecy antecedent
  Term new_prop =
      fts.make_term(Implies, fts.make_term(Equal, proph_var, new_target), prop);

  UnorderedTermSet free_vars;
  get_free_symbolic_consts(new_prop, free_vars);

  // Expecting two new history variable and one new prophecy variable
  EXPECT_EQ(fts.statevars().size() - 3, num_statevars_orig);

  // x should still be in the property
  EXPECT_TRUE(free_vars.find(x) != free_vars.end());

  // but now the prophecy variable should be also
  EXPECT_TRUE(free_vars.find(proph_var) != free_vars.end());
}

TEST_P(ModifierUnitTests, ImplicitPredicateAbstractor)
{
  RelationalTransitionSystem rts(s);
  counter_system(rts, rts.make_term(10, bvsort));
  Term x = rts.named_terms().at("x");

  RelationalTransitionSystem abs_rts(rts.solver());
  Unroller un(abs_rts);
  ImplicitPredicateAbstractor ia(rts, abs_rts, un);

  ia.do_abstraction();

  // check if c <= 10 is inductive on the concrete system
  Term x_le_10 = rts.make_term(BVUle, x, rts.make_term(10, bvsort));
  s->push();
  s->assert_formula(x_le_10);
  s->assert_formula(rts.trans());
  s->assert_formula(s->make_term(Not, rts.next(x_le_10)));
  Result r = s->check_sat();
  s->pop();
  EXPECT_TRUE(r.is_unsat());  // expecting it to be inductive

  // check if c <= 10 is inductive on the abstract system
  s->push();
  s->assert_formula(x_le_10);
  s->assert_formula(abs_rts.trans());
  s->assert_formula(s->make_term(Not, abs_rts.next(x_le_10)));
  r = s->check_sat();
  s->pop();
  EXPECT_TRUE(r.is_sat());  // expecting check to fail

  // add it as a predicate
  Term ref = ia.predicate_refinement(x_le_10);
  abs_rts.constrain_trans(ref);

  // check if c <= 10 is inductive on the refined abstract system
  s->push();
  s->assert_formula(x_le_10);
  s->assert_formula(abs_rts.trans());
  s->assert_formula(s->make_term(Not, abs_rts.next(x_le_10)));
  r = s->check_sat();
  s->pop();
  EXPECT_TRUE(r.is_unsat());  // expecting it to be inductive now
}

TEST_P(ModifierUnitTests, LivenessToSafetyTranslator)
{
  // Create TS with one counter state x mod 4.
  FunctionalTransitionSystem fts(s);
  Term max_val = fts.make_term(3, bvsort);
  counter_system(fts, max_val);
  EXPECT_EQ(fts.inputvars().size(), 0);
  EXPECT_EQ(fts.statevars().size(), 1);

  // Create justice property x = -1.
  Term x = fts.named_terms().at("x");
  Term minus_one = fts.make_term(-1, bvsort);
  Term justice_cond = fts.make_term(smt::Equal, x, minus_one);

  // Perform liveness to safety translation.
  LivenessToSafetyTranslator l2s;
  Term prop = l2s.translate(fts, { justice_cond });

  // Expect one new boolean input "save".
  EXPECT_EQ(fts.inputvars().size(), 1);
  Term save;
  for (auto inputvar : fts.inputvars()) {
    save = inputvar;
  }
  EXPECT_EQ(save->get_sort(), boolsort);

  // Expect three new states: "saved" (bool), "justice" (bool), "loop" (bv).
  EXPECT_EQ(fts.statevars().size(), 4);
  Term saved;
  Term justice;
  Term loop;
  for (auto statevar : fts.statevars()) {
    if (statevar == x) continue;
    if (statevar->get_sort() == boolsort) {
      auto state_name = statevar->to_string();
      if (state_name.find("saved") != state_name.npos) {
        saved = statevar;
      } else {
        justice = statevar;
      }
    } else if (statevar->get_sort() == bvsort) {
      loop = statevar;
    } else {
      FAIL() << "unexpected sort for state var";
    }
  }

  // saved' = save | saved
  Term saved_expected = fts.make_term(smt::Or, saved, save);
  Term saved_actual = fts.state_updates().at(saved);
  Term saved_correct = fts.make_term(smt::Equal, saved_actual, saved_expected);
  s->push();
  s->assert_formula(fts.make_term(smt::Not, saved_correct));
  auto result = s->check_sat();
  s->pop();
  EXPECT_TRUE(result.is_unsat());

  // loop' = (!saved & save) ? x : loop
  Term loop_expected = fts.make_term(
      smt::Ite,
      { fts.make_term(smt::And, fts.make_term(smt::Not, saved), save),
        x,
        loop });
  Term loop_actual = fts.state_updates().at(loop);
  Term loop_correct = fts.make_term(smt::Equal, loop_actual, loop_expected);
  s->push();
  s->assert_formula(fts.make_term(smt::Not, loop_correct));
  result = s->check_sat();
  s->pop();
  EXPECT_TRUE(result.is_unsat());

  // Safety property: !((saved & x = loop) & justice)
  Term looped =
      fts.make_term(smt::And, saved, fts.make_term(smt::Equal, x, loop));
  Term bad = fts.make_term(smt::And, looped, justice);
  Term safety = fts.make_term(smt::Not, bad);
  Term prop_correct = fts.make_term(smt::Equal, prop, safety);
  s->push();
  s->assert_formula(fts.make_term(smt::Not, prop_correct));
  result = s->check_sat();
  s->pop();
  EXPECT_TRUE(result.is_unsat());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedModifierUnitTests,
                         ModifierUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
