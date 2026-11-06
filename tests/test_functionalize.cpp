#include <vector>

#include "core/fts.h"
#include "core/prop.h"
#include "core/proverresult.h"
#include "core/rts.h"
#include "core/ts.h"
#include "engines/bmc.h"
#include "gtest/gtest.h"
#include "modifiers/functionalize.h"
#include "smt-switch/smt.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/str_util.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

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
  Bmc bmc(SafetyProperty(s, prop), fts, prover_solver(BMC));
  ASSERT_EQ(bmc.check_until(10), FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(bmc.witness(cex));
  ASSERT_EQ(cex.size(), 4);
  for (size_t k = 0; k < cex.size(); ++k) {
    EXPECT_EQ(cex[k].at(x), bv(static_cast<int64_t>(k)));
  }
}

INSTANTIATE_TEST_SUITE_P(ParameterizedFunctionalizeTests,
                         FunctionalizeTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
