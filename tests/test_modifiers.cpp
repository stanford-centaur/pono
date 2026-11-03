#include <tuple>
#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/rts.h"
#include "engines/bmc.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "modifiers/array_abstractor.h"
#include "modifiers/history_modifier.h"
#include "modifiers/implicit_predicate_abstractor.h"
#include "modifiers/liveness_to_safety_translator.h"
#include "modifiers/mod_ts_prop.h"
#include "modifiers/prophecy_modifier.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/ts_analysis.h"

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
    s = create_solver(GetParam(), GetParam() == BTOR);
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

TEST_P(ModifierUnitTests, HistoryModifierHierarchicalName)
{
  // Regression test: get_hist() builds each history variable's name by
  // concatenating a prefix/suffix onto Term::to_string(), which needs
  // desanitizing for names that require SMT-LIB `|...|` quoting (e.g.
  // hierarchical names like SystemVerilog produces, such as "mod.sig[3]").
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("mod.sig[3]", bvsort);
  fts.assign_next(x, x);

  HistoryModifier hm(fts);
  ASSERT_NO_THROW(hm.get_hist(x, 3));
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

TEST_P(ModifierUnitTests, ProphecyModifierHierarchicalName)
{
  // Regression test: get_proph() builds the prophecy variable's name by
  // concatenating a prefix/suffix onto Term::to_string(), which needs
  // desanitizing for names that require SMT-LIB `|...|` quoting. It also
  // composes HistoryModifier::get_hist() internally, so this exercises
  // both fixes together.
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("mod.sig[3]", bvsort);
  fts.assign_next(x, x);

  ProphecyModifier pm(fts);
  ASSERT_NO_THROW(pm.get_proph(x, 2));
}

TEST_P(ModifierUnitTests, ImplicitPredicateAbstractorHierarchicalName)
{
  // Regression test: do_abstraction() builds the abstracted next-var's name
  // by concatenating "^" onto Term::to_string(), which needs desanitizing
  // for names that require SMT-LIB `|...|` quoting.
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("mod.sig[3]", bvsort);
  rts.assign_next(x, x);
  rts.set_init(rts.make_term(Equal, x, rts.make_term(0, bvsort)));

  RelationalTransitionSystem abs_rts(rts.solver());
  Unroller un(abs_rts);
  ImplicitPredicateAbstractor ia(rts, abs_rts, un);
  ASSERT_NO_THROW(ia.do_abstraction());
}

TEST_P(ModifierUnitTests, ArrayAbstractorHierarchicalName)
{
  if (s->get_solver_enum() == BTOR) {
    GTEST_SKIP() << "Boolector does not support abstract sorts";
  }
  if (!solver_has_attribute(GetParam(), smt::CONSTARR)) {
    GTEST_SKIP() << GetParam() << " does not support constant arrays";
  }
  // Regression test: ArrayAbstractor builds names for the abstracted
  // constant-array, state, and input variables by concatenating a prefix
  // onto Term::to_string(), which needs desanitizing for names that require
  // SMT-LIB `|...|` quoting.
  RelationalTransitionSystem conc_rts(s);
  Term arr = conc_rts.make_statevar("mod.arr[3]", arrsort);
  Term arr_in = conc_rts.make_inputvar("mod.in_arr[1]", arrsort);
  Term constarr0 = conc_rts.make_term(conc_rts.make_term(0, bvsort), arrsort);
  conc_rts.set_init(conc_rts.make_term(Equal, arr, constarr0));
  conc_rts.assign_next(arr, arr);

  RelationalTransitionSystem abs_rts(s);
  ArrayAbstractor aa(conc_rts, abs_rts);
  ASSERT_NO_THROW(aa.do_abstraction());
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

TEST_P(ModifierUnitTests, LivenessToSafetyTranslatorHierarchicalName)
{
  // Regression test: a state variable whose name requires SMT-LIB `|...|`
  // quoting (e.g. a hierarchical name like SystemVerilog produces, such as
  // "mod.sig[3]") must not corrupt the generated "loop" copy's name.
  // LivenessToSafetyTranslator used to build the new name via
  // `statevar->to_string() + suffix`; to_string() already returns the
  // quoted form for such names, so appending more text after it produced a
  // malformed name (text trailing a closing quote), which some solvers
  // (e.g. bzla) reject with "invalid symbol ... not SMT-LIB compliant".
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("mod.sig[3]", bvsort);
  fts.assign_next(x, x);
  Term minus_one = fts.make_term(-1, bvsort);
  Term justice_cond = fts.make_term(smt::Equal, x, minus_one);

  LivenessToSafetyTranslator l2s;
  ASSERT_NO_THROW(l2s.translate(fts, { justice_cond }));

  Term loop;
  for (auto statevar : fts.statevars()) {
    if (statevar != x && statevar->get_sort() == bvsort) {
      loop = statevar;
    }
  }
  ASSERT_TRUE(loop);

  // The printed name must be a single well-formed SMT-LIB token: after
  // stripping one optional outer `|...|` quote pair, there must be no
  // leftover '|' characters.
  string loop_name = loop->to_string();
  if (loop_name.size() > 2 && loop_name.front() == '|'
      && loop_name.back() == '|') {
    loop_name = loop_name.substr(1, loop_name.size() - 2);
  }
  EXPECT_EQ(loop_name.find('|'), string::npos);
}

TEST_P(ModifierUnitTests, EliminateFloorAndDivWithoutIntegers)
{
  // a system w/o integers has nothing to replace and is left as it is
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  rts.set_init(rts.make_term(Equal, x, rts.make_term(0, bvsort)));
  rts.assign_next(x, rts.make_term(BVAdd, x, rts.make_term(1, bvsort)));
  Term prop = rts.make_term(BVUlt, x, rts.make_term(10, bvsort));
  Term init = rts.init();
  Term trans = rts.trans();
  eliminate_floor_and_div(rts, prop);
  EXPECT_EQ(rts.init(), init);
  EXPECT_EQ(rts.trans(), trans);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedModifierUnitTests,
                         ModifierUnitTests,
                         testing::ValuesIn(available_solver_enums()));

// the base of the suites below, each instantiated with the solvers to test,
// such as cvc5 and MathSAT: a suite builds a system with div, mod or floor on
// one of them and checks that BMC finds the same counterexample in it, or
// k-induction proves it at the same bound, before and after the elimination
class EliminateFloorAndDivTest : public ::testing::Test
{
 protected:
  // whether a term contains integer div or mod, or floor
  static bool has_div_mod(const Term & t)
  {
    UnorderedOpSet ops;
    get_ops(t, ops);
    return ops.count(Op(IntDiv)) || ops.count(Op(Mod)) || ops.count(Op(To_Int));
  }

  // starts rts, empty, on a new solver of kind se
  void start_system(SolverEnum se)
  {
    solver_enum = se;
    s = create_solver(se);
    intsort = s->make_sort(INT);
    rts = RelationalTransitionSystem(s);
  }

  // makes rts the system of DivModTests: x counts up from 0, and y is the
  // previous x mod 3
  void build_div_mod_system()
  {
    x = rts.make_statevar("x", intsort);
    y = rts.make_statevar("y", intsort);
    three = s->make_term(3, intsort);
    Term zero = s->make_term(0, intsort);
    rts.set_init(s->make_term(
        And, s->make_term(Equal, x, zero), s->make_term(Equal, y, zero)));
    rts.assign_next(x, s->make_term(Plus, x, s->make_term(1, intsort)));
    rts.assign_next(y, s->make_term(Mod, x, three));
  }

  // makes rts the system of FloorTests: r counts up from 0 in steps of 1/2,
  // and f is the previous floor(r)
  void build_floor_system()
  {
    Sort realsort = s->make_sort(REAL);
    r = rts.make_statevar("r", realsort);
    f = rts.make_statevar("f", intsort);
    rts.set_init(
        s->make_term(And,
                     s->make_term(Equal, r, s->make_term(0, realsort)),
                     s->make_term(Equal, f, s->make_term(0, intsort))));
    Term half =
        s->make_term(Div, s->make_term(1, realsort), s->make_term(2, realsort));
    rts.assign_next(r, s->make_term(Plus, r, half));
    rts.assign_next(f, s->make_term(To_Int, r));
  }

  // a fresh solver of the kind under test, for a prover: each prover declares
  // its own symbols, which a second one on the same solver would clash with
  SmtSolver fresh_solver() const { return create_solver(solver_enum); }

  // a copy of rts and prop with div, mod and floor eliminated
  pair<RelationalTransitionSystem, Term> eliminated(
      Term prop, bool allow_new_inputvars = true) const
  {
    RelationalTransitionSystem elim = rts;
    eliminate_floor_and_div(elim, prop, allow_new_inputvars);
    return { elim, prop };
  }

  // the counterexample that BMC finds for prop in ts within 10 steps, as the
  // values of the variables of ts in each step, or none if there is none
  vector<UnorderedTermMap> bmc_cex(const TransitionSystem & ts,
                                   const Term & prop) const
  {
    Bmc bmc({ s, prop }, ts, fresh_solver());
    vector<UnorderedTermMap> cex;
    if (bmc.check_until(10) == ProverResult::FALSE) {
      EXPECT_TRUE(bmc.witness(cex));
    }
    return cex;
  }

  // the least bound at which k-induction proves prop in ts, or -1 if there is
  // none up to 10
  int kind_proof_bound(const TransitionSystem & ts, const Term & prop) const
  {
    for (int k = 0; k <= 10; k++) {
      KInduction kind({ s, prop }, ts, fresh_solver());
      if (kind.check_until(k) == ProverResult::TRUE) {
        return k;
      }
    }
    return -1;
  }

  // checks that BMC finds a counterexample to prop in the original rts, and
  // the same one after elimination. The counterexample must be unique, as it
  // is in the systems built for these tests: init fixes each variable and its
  // update is a function of the state.
  void expect_same_cex(const Term & prop, bool allow_new_inputvars = true) const
  {
    auto [elim, elim_prop] = eliminated(prop, allow_new_inputvars);
    vector<UnorderedTermMap> cex = bmc_cex(rts, prop);
    ASSERT_FALSE(cex.empty());
    vector<UnorderedTermMap> elim_cex = bmc_cex(elim, elim_prop);
    ASSERT_EQ(elim_cex.size(), cex.size());
    for (size_t i = 0; i < cex.size(); i++) {
      for (const auto & v : rts.statevars()) {
        EXPECT_EQ(elim_cex[i].at(v)->to_string(), cex[i].at(v)->to_string())
            << v << " in step " << i;
      }
    }
  }

  // checks that k-induction proves prop in the original rts, and at the same
  // bound after elimination
  void expect_same_kind_proof_bound(const Term & prop) const
  {
    auto [elim, elim_prop] = eliminated(prop);
    int k = kind_proof_bound(rts, prop);
    ASSERT_GE(k, 0);
    EXPECT_EQ(kind_proof_bound(elim, elim_prop), k);
  }

  SolverEnum solver_enum;
  SmtSolver s;
  Sort intsort;
  RelationalTransitionSystem rts;
  Term x, y, three;  // in the system of DivModTests
  Term r, f;         // in the system of FloorTests
};

class DivModTests : public EliminateFloorAndDivTest,
                    public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    start_system(GetParam());
    build_div_mod_system();
  }
};

TEST_P(DivModTests, NoneLeft)
{
  // also a div over next-state variables, and a mod in the property
  Term next_x = rts.next(x);
  rts.constrain_trans(
      s->make_term(Le, s->make_term(IntDiv, next_x, three), next_x));
  Term prop = s->make_term(Lt, s->make_term(Mod, x, three), three);
  ASSERT_TRUE(has_div_mod(rts.trans()));
  ASSERT_TRUE(has_div_mod(prop));

  eliminate_floor_and_div(rts, prop);
  EXPECT_FALSE(has_div_mod(rts.init()));
  EXPECT_FALSE(has_div_mod(rts.trans()));
  EXPECT_FALSE(has_div_mod(prop));
  // x, y and one new variable q for (div x 3): (mod x 3) becomes x - 3*q, and
  // (div x' 3) becomes q', the next-state copy of q
  EXPECT_EQ(rts.statevars().size(), 3);
}

TEST_P(DivModTests, CurrentAndNextState)
{
  // z' is (div (+ x y') 3), over both current- and next-state variables, so
  // it gets a variable of its own rather than the next-state copy of the one
  // for (div (+ x y) 3). z is 1 first in step 3.
  Term z = rts.make_statevar("z", intsort);
  rts.constrain_init(s->make_term(Equal, z, s->make_term(0, intsort)));
  rts.constrain_trans(s->make_term(
      Equal,
      rts.next(z),
      s->make_term(IntDiv, s->make_term(Plus, x, rts.next(y)), three)));
  Term prop = s->make_term(Lt, z, s->make_term(1, intsort));

  auto [elim, elim_prop] = eliminated(prop);
  EXPECT_FALSE(has_div_mod(elim.trans()));
  // x, y, z and a new state variable for (div x 3), and an input variable for
  // (div (+ x y') 3), as its bounds constrain only the transition
  EXPECT_EQ(elim.statevars().size(), 4);
  EXPECT_EQ(elim.inputvars().size(), 1);
  expect_same_cex(prop);

  // a state variable for (div (+ x y') 3) instead, if no new input variables
  // are allowed, as with --promote-inputvars
  auto [elim_sv, elim_sv_prop] = eliminated(prop, false);
  EXPECT_EQ(elim_sv.statevars().size(), 5);
  EXPECT_TRUE(elim_sv.inputvars().empty());
  expect_same_cex(prop, false);
}

TEST_P(DivModTests, Safe)
{
  expect_same_kind_proof_bound(s->make_term(Lt, y, three));
}

TEST_P(DivModTests, Unsafe)
{
  // y is 2 first in step 3
  expect_same_cex(s->make_term(Lt, y, s->make_term(2, intsort)));
}

TEST_P(DivModTests, NegatedDivisor)
{
  // a negated literal, which cvc5 and Z3 keep as (- 3)
  Term prop = s->make_term(
      Lt, s->make_term(Mod, x, s->make_term(Negate, three)), three);
  eliminate_floor_and_div(rts, prop);
  EXPECT_FALSE(has_div_mod(prop));
}

TEST_P(DivModTests, ZeroDivisor)
{
  // (- 3 3) and (- 0) have no symbols but are 0, so div and mod by them must
  // be kept: c*q <= t < c*q + |c| has no solution for c = 0.
  Term zero = s->make_term(0, intsort);
  Term prop = s->make_term(
      And,
      s->make_term(
          Le, s->make_term(IntDiv, x, s->make_term(Minus, three, three)), x),
      s->make_term(Le, s->make_term(Mod, x, s->make_term(Negate, zero)), x));
  UnorderedOpSet before, after;
  get_ops(prop, before);
  eliminate_floor_and_div(rts, prop);
  get_ops(prop, after);
  EXPECT_EQ(after.count(Op(IntDiv)), before.count(Op(IntDiv)));
  EXPECT_EQ(after.count(Op(Mod)), before.count(Op(Mod)));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedDivModTests,
    DivModTests,
    testing::ValuesIn(filter_solver_enums({ THEORY_INT })));

class FloorTests : public EliminateFloorAndDivTest,
                   public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    start_system(GetParam());
    build_floor_system();
  }
};

TEST_P(FloorTests, NoneLeft)
{
  Term prop = s->make_term(Le, s->make_term(To_Real, f), r);
  ASSERT_TRUE(has_div_mod(rts.trans()));
  auto [elim, elim_prop] = eliminated(prop);
  EXPECT_FALSE(has_div_mod(elim.init()));
  EXPECT_FALSE(has_div_mod(elim.trans()));
}

TEST_P(FloorTests, Safe)
{
  expect_same_kind_proof_bound(s->make_term(Le, s->make_term(To_Real, f), r));
}

TEST_P(FloorTests, Unsafe)
{
  // f is 2 first in step 5, the floor of r = 2 in step 4
  expect_same_cex(s->make_term(Lt, f, s->make_term(2, intsort)));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedFloorTests,
    FloorTests,
    testing::ValuesIn(filter_solver_enums({ THEORY_INT, THEORY_REAL })));

// the systems of DivModTests and FloorTests, proved by interpolation with an
// interpolator under test as well: MathSAT's interpolator fails on floor,
// which smt-switch uses to encode div and mod for MathSAT
class DivModInterpTests
    : public EliminateFloorAndDivTest,
      public ::testing::WithParamInterface<tuple<SolverEnum, SolverEnum>>
{
 protected:
  void SetUp() override { start_system(get<0>(GetParam())); }

  // checks that interpolation proves prop after elimination, with an
  // invariant that holds in the original rts as well once the one new
  // variable in it is replaced by the term the variable stands for
  void expect_valid_imc_invariant(const Term & prop,
                                  const Term & replaced) const
  {
    auto [elim, elim_prop] = eliminated(prop);
    UnorderedTermMap var_to_term;
    for (const auto & v : elim.statevars()) {
      if (!rts.is_curr_var(v)) {
        var_to_term[v] = replaced;
      }
    }
    ASSERT_EQ(var_to_term.size(), 1);

    PonoOptions opts;
    opts.smt_interpolator_ = get<1>(GetParam());
    InterpolantMC imc({ s, elim_prop }, elim, fresh_solver(), opts);
    ASSERT_EQ(imc.check_until(10), ProverResult::TRUE);
    Term invar = s->substitute(imc.invar(), var_to_term);
    EXPECT_TRUE(check_invar(rts, prop, invar));
  }
};

TEST_P(DivModInterpTests, DivMod)
{
  build_div_mod_system();
  expect_valid_imc_invariant(s->make_term(Lt, y, three),
                             s->make_term(IntDiv, x, three));
}

TEST_P(DivModInterpTests, Floor)
{
  if (get<1>(GetParam()) == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5's interpolator gives up on this system";
  }
  build_floor_system();
  expect_valid_imc_invariant(s->make_term(Le, s->make_term(To_Real, f), r),
                             s->make_term(To_Int, r));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedDivModInterpTests,
    DivModInterpTests,
    testing::Combine(
        testing::ValuesIn(filter_solver_enums({ THEORY_INT, THEORY_REAL })),
        testing::ValuesIn(filter_interpolator_enums({ THEORY_INT }))));

}  // namespace pono_tests
