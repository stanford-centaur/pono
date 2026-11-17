#include <vector>

#include "core/fts.h"
#include "core/prop.h"
#include "core/proverresult.h"
#include "core/rts.h"
#include "core/ts.h"
#include "engines/bmc.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "modifiers/functionalize.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/str_util.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

/** @return the conjuncts of a term, regardless of order and repetition */
UnorderedTermSet conjunct_set(const Term & term)
{
  TermVec conjuncts;
  conjunctive_partition(term, conjuncts, true);
  return UnorderedTermSet(conjuncts.begin(), conjuncts.end());
}

class FunctionalizeTests : public ::testing::Test,
                           public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam(), GetParam() == BTOR);
    bvsort = s->make_sort(BV, 8);
  }

  Term bv(int64_t value) { return s->make_term(value, bvsort); }

  /** @return a solver of its own for a prover, whose unroller names the
   *  timed copies of the variables after them */
  SmtSolver prover_solver(Engine e)
  {
    return create_solver_for(GetParam(), e, GetParam() == BTOR);
  }

  /** Checks a property of a system with bounded model checking.
   *  @return the length of the counterexample, or -1 if there is none
   */
  int64_t get_cex_len(const TransitionSystem & ts, const Term & prop)
  {
    Bmc bmc(SafetyProperty(s, prop), ts, prover_solver(BMC));
    if (bmc.check_until(10) != FALSE) {
      return -1;
    }
    return static_cast<int64_t>(bmc.witness_length());
  }

  SmtSolver s;
  Sort bvsort;
};

TEST_P(FunctionalizeTests, KeepsFunctionalSystem)
{
  FunctionalTransitionSystem fts(s);
  counter_system(fts, bv(10));
  EXPECT_EQ(functionalize(fts), fts);
}

// The variable x_in that chooses the next value of x is a new input, or,
// without new inputs, a state variable without an update, which takes any
// value in each step as well.
TEST_P(FunctionalizeTests, ChoosesNextValuesByInputsOrStateVariables)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term idle = rts.make_statevar("idle", bvsort);
  Term in = rts.make_inputvar("in", bvsort);
  rts.set_init(rts.make_term(Equal, x, bv(0)));
  rts.constrain_trans(
      rts.make_term(BVUgt, rts.next(x), rts.make_term(BVAdd, x, in)));
  rts.constrain_trans(rts.make_term(BVUlt, rts.next(x), bv(10)));
  rts.add_constraint(rts.make_term(BVUlt, in, bv(5)));

  for (bool allow_new_inputvars : { true, false }) {
    SCOPED_TRACE(allow_new_inputvars ? "new inputs" : "no new inputs");
    TransitionSystem fts = functionalize(rts, allow_new_inputvars);
    ASSERT_TRUE(fts.is_functional());
    Term x_in = fts.state_updates().at(x);
    EXPECT_TRUE(is_generated_name(name_desanitize(x_in->to_string())));
    // idle takes any value in the next state, as before
    UnorderedTermSet inputs = { in };
    UnorderedTermSet no_update = { idle };
    (allow_new_inputvars ? inputs : no_update).insert(x_in);
    EXPECT_EQ(fts.inputvars(), inputs);
    EXPECT_EQ(fts.statevars_with_no_update(), no_update);
    // the updates of x and valid, and the constraint and valid
    EXPECT_EQ(fts.state_updates().size(), 2);
    EXPECT_EQ(fts.constraints().size(), 2);
    EXPECT_EQ(fts.lookup("x"), x);
    for (int64_t value : { 5, 9, 10 }) {
      Term prop = rts.make_term(Distinct, x, bv(value));
      EXPECT_EQ(get_cex_len(fts, prop), get_cex_len(rts, prop))
          << "reaching " << value;
    }
  }
}

// A counterexample keeps its length and values, even if its last state has
// no successor, as valid checks each step in the state it leads to.
TEST_P(FunctionalizeTests, PreservesCex)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  rts.set_init(rts.make_term(Equal, x, bv(0)));
  rts.constrain_trans(rts.make_term(BVUlt, x, bv(3)));
  rts.constrain_trans(
      rts.make_term(Equal, rts.next(x), rts.make_term(BVAdd, x, bv(1))));
  Term prop = rts.make_term(Distinct, x, bv(3));
  ASSERT_EQ(get_cex_len(rts, prop), 3);

  TransitionSystem fts = functionalize(rts);
  // the equality gives x its update, and valid checks the guard
  EXPECT_TRUE(fts.inputvars().empty());
  EXPECT_EQ(fts.statevars().size(), 2);
  Bmc bmc(SafetyProperty(s, prop), fts, prover_solver(BMC));
  ASSERT_EQ(bmc.check_until(10), FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(bmc.witness(cex));
  ASSERT_EQ(cex.size(), 4);
  for (size_t k = 0; k < cex.size(); ++k) {
    EXPECT_EQ(cex[k].at(x), bv(static_cast<int64_t>(k)));
  }
}

TEST_P(FunctionalizeTests, TakesUpdatesFromEqualities)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term y = rts.make_statevar("y", bvsort);
  Term z = rts.make_statevar("z", bvsort);
  Term in = rts.make_inputvar("in", bvsort);
  for (const Term & v : { x, y, z }) {
    rts.constrain_init(rts.make_term(Equal, v, bv(0)));
  }
  // z follows the next values of x and y, which come later
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(z), rts.make_term(BVAdd, rts.next(x), rts.next(y))));
  rts.constrain_trans(
      rts.make_term(Equal, rts.next(x), rts.make_term(BVAdd, x, in)));
  rts.constrain_trans(
      rts.make_term(Equal, rts.make_term(BVAdd, y, bv(1)), rts.next(y)));

  TransitionSystem fts = functionalize(rts);
  // As every conjunct defines an update, no input chooses a next value, and
  // there is no valid.
  EXPECT_EQ(fts.inputvars(), UnorderedTermSet({ in }));
  EXPECT_EQ(fts.statevars().size(), 3);
  EXPECT_TRUE(fts.statevars_with_no_update().empty());
  EXPECT_TRUE(fts.constraints().empty());
  EXPECT_EQ(fts.state_updates().at(x), rts.make_term(BVAdd, x, in));
  EXPECT_EQ(fts.state_updates().at(y), rts.make_term(BVAdd, y, bv(1)));

  Term sum = rts.make_term(Equal, z, rts.make_term(BVAdd, x, y));
  KInduction kind(SafetyProperty(s, sum), fts, prover_solver(KIND));
  EXPECT_EQ(kind.check_until(10), TRUE);
  Term prop = rts.make_term(Distinct, z, bv(5));
  EXPECT_EQ(get_cex_len(fts, prop), get_cex_len(rts, prop));
}

// A functional system seen as a relational one comes back as it was: the
// updates come from the equalities, the constraints stay constraints, and no
// auxiliary variable is needed. Only the order of the conjuncts may differ,
// and the initial states get the constraints once more.
TEST_P(FunctionalizeTests, RoundTripGetsSameFunctionalSystem)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  Term y = fts.make_statevar("y", bvsort);
  Term in = fts.make_inputvar("in", bvsort);
  fts.constrain_init(fts.make_term(Equal, x, bv(0)));
  fts.constrain_init(fts.make_term(Equal, y, bv(0)));
  fts.assign_next(x, fts.make_term(BVAdd, x, in));
  fts.assign_next(y, fts.make_term(Ite, fts.make_term(BVUlt, x, y), y, x));
  fts.add_constraint(fts.make_term(BVUlt, in, bv(5)));
  fts.add_constraint(fts.make_term(Distinct, x, bv(200)));
  fts.name_term("sum", fts.make_term(BVAdd, x, y));

  TransitionSystem round_trip = functionalize(RelationalTransitionSystem(fts));
  // not round_trip == fts, which compares init and trans as terms, but
  // functionalize may reorder and repeat their conjuncts
  EXPECT_TRUE(round_trip.is_functional());
  EXPECT_EQ(round_trip.statevars(), fts.statevars());
  EXPECT_EQ(round_trip.inputvars(), fts.inputvars());
  EXPECT_EQ(round_trip.state_updates(), fts.state_updates());
  EXPECT_EQ(round_trip.constraints(), fts.constraints());
  EXPECT_EQ(round_trip.named_terms(), fts.named_terms());
  EXPECT_EQ(conjunct_set(round_trip.init()), conjunct_set(fts.init()));
  EXPECT_EQ(conjunct_set(round_trip.trans()), conjunct_set(fts.trans()));
}

// Equalities that define next-state variables by each other depend on each
// other, until an input becomes the update of one of the variables. A chain
// of equalities that leads into the cycle needs no input of its own.
TEST_P(FunctionalizeTests, CutsCycleOfEqualities)
{
  RelationalTransitionSystem rts(s);
  Term u = rts.make_statevar("u", bvsort);
  Term x = rts.make_statevar("x", bvsort);
  Term y = rts.make_statevar("y", bvsort);
  Term z = rts.make_statevar("z", bvsort);
  rts.constrain_init(rts.make_term(Equal, u, bv(1)));
  for (const Term & v : { x, y, z }) {
    rts.constrain_init(rts.make_term(Equal, v, bv(0)));
  }
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(u), rts.make_term(BVAdd, rts.next(x), bv(1))));
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(x), rts.make_term(BVAdd, rts.next(y), rts.next(z))));
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(y), rts.make_term(BVSub, rts.next(x), rts.next(z))));
  rts.constrain_trans(
      rts.make_term(Equal, rts.next(z), rts.make_term(BVAdd, z, bv(1))));

  TransitionSystem fts = functionalize(rts);
  EXPECT_EQ(fts.state_updates().at(z), rts.make_term(BVAdd, z, bv(1)));
  // an input chooses the next value of x or y, which defines the other one,
  // and valid checks the equality for the first
  EXPECT_EQ(fts.inputvars().size(), 1);
  EXPECT_NE(fts.is_input_var(fts.state_updates().at(x)),
            fts.is_input_var(fts.state_updates().at(y)));
  EXPECT_FALSE(fts.is_input_var(fts.state_updates().at(u)));
  EXPECT_EQ(fts.constraints().size(), 1);

  Term follows =
      rts.make_term(And,
                    rts.make_term(Equal, x, rts.make_term(BVAdd, y, z)),
                    rts.make_term(Equal, u, rts.make_term(BVAdd, x, bv(1))));
  KInduction kind(SafetyProperty(s, follows), fts, prover_solver(KIND));
  EXPECT_EQ(kind.check_until(10), TRUE);
  for (int64_t value : { 7, 200 }) {
    Term prop = rts.make_term(Distinct, x, bv(value));
    EXPECT_EQ(get_cex_len(fts, prop), get_cex_len(rts, prop))
        << "reaching " << value;
  }

  // without new inputs, a state variable without an update cuts the cycle
  TransitionSystem no_inputs = functionalize(rts, false);
  EXPECT_TRUE(no_inputs.inputvars().empty());
  EXPECT_EQ(no_inputs.statevars_with_no_update().size(), 1);
  Term prop = rts.make_term(Distinct, x, bv(7));
  EXPECT_EQ(get_cex_len(no_inputs, prop), get_cex_len(rts, prop));
}

// An equality that depends on a next-state variable no equality defines is
// resolved by the input that the variable gets anyway.
TEST_P(FunctionalizeTests, CutsAtVariableNoEqualityDefines)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term y = rts.make_statevar("y", bvsort);
  rts.constrain_init(rts.make_term(Equal, x, bv(1)));
  rts.constrain_init(rts.make_term(Equal, y, bv(0)));
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(x), rts.make_term(BVAdd, rts.next(y), bv(1))));
  rts.constrain_trans(rts.make_term(BVUlt, rts.next(y), bv(5)));

  TransitionSystem fts = functionalize(rts);
  EXPECT_EQ(fts.inputvars().size(), 1);
  EXPECT_TRUE(fts.is_input_var(fts.state_updates().at(y)));
  EXPECT_FALSE(fts.is_input_var(fts.state_updates().at(x)));
  // valid checks the bound on y
  EXPECT_EQ(fts.constraints().size(), 1);

  Term follows = rts.make_term(Equal, x, rts.make_term(BVAdd, y, bv(1)));
  KInduction kind(SafetyProperty(s, follows), fts, prover_solver(KIND));
  EXPECT_EQ(kind.check_until(10), TRUE);
  for (int64_t value : { 5, 6 }) {
    Term prop = rts.make_term(Distinct, x, bv(value));
    EXPECT_EQ(get_cex_len(fts, prop), get_cex_len(rts, prop))
        << "reaching " << value;
  }
}

// An equality that defines a variable of a cycle by current-state variables
// lets the other equalities of the cycle define updates.
TEST_P(FunctionalizeTests, TakesUpdatesFromCycleAnotherEqualityBreaks)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term y = rts.make_statevar("y", bvsort);
  Term z = rts.make_statevar("z", bvsort);
  rts.constrain_init(rts.make_term(Equal, x, bv(0)));
  rts.constrain_init(rts.make_term(Equal, y, bv(1)));
  rts.constrain_init(rts.make_term(Equal, z, bv(0)));
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(y), rts.make_term(BVAdd, rts.next(x), bv(1))));
  rts.constrain_trans(rts.make_term(
      Equal, rts.next(x), rts.make_term(BVSub, rts.next(y), bv(1))));
  rts.constrain_trans(
      rts.make_term(Equal, rts.next(x), rts.make_term(BVAdd, z, bv(1))));
  rts.constrain_trans(
      rts.make_term(Equal, rts.next(z), rts.make_term(BVAdd, z, bv(1))));

  TransitionSystem fts = functionalize(rts);
  EXPECT_TRUE(fts.inputvars().empty());
  EXPECT_TRUE(fts.statevars_with_no_update().empty());
  EXPECT_EQ(fts.state_updates().at(x), rts.make_term(BVAdd, z, bv(1)));
  // valid checks the other equality for x, x' = y' - 1
  EXPECT_EQ(fts.constraints().size(), 1);

  Term follows = rts.make_term(Equal, y, rts.make_term(BVAdd, x, bv(1)));
  KInduction kind(SafetyProperty(s, follows), fts, prover_solver(KIND));
  EXPECT_EQ(kind.check_until(10), TRUE);
  Term prop = rts.make_term(Distinct, y, bv(3));
  EXPECT_EQ(get_cex_len(fts, prop), get_cex_len(rts, prop));
}

INSTANTIATE_TEST_SUITE_P(ParameterizedFunctionalizeTests,
                         FunctionalizeTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
