#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/rts.h"
#include "engines/bmc.h"
#include "modifiers/control_signals.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class ControlUnitTests : public ::testing::Test,
                         public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s->set_opt("produce-models", "true");
    s->set_opt("incremental", "true");
    boolsort = s->make_sort(BOOL);
    bvsort1 = s->make_sort(BV, 1);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort1, bvsort8;
};

TEST_P(ControlUnitTests, SimpleReset)
{
  FunctionalTransitionSystem fts(s);
  Term rst = fts.make_inputvar("rst", boolsort);
  Term x = fts.make_statevar("x", bvsort8);
  // x' = rst ? 0 : (x < 10) ? x + 1 : x
  Term x_update =
      fts.make_term(Ite,
                    fts.make_term(BVUlt, x, fts.make_term(10, bvsort8)),
                    fts.make_term(BVAdd, x, fts.make_term(1, bvsort8)),
                    x);
  fts.assign_next(x,
                  fts.make_term(Ite, rst, fts.make_term(0, bvsort8), x_update));
  Term p_false_term = fts.make_term(BVUle, x, fts.make_term(10, bvsort8));
  Property p_false(fts, p_false_term);
  // use a new context so unroller doesn't clash on unrolled symbols
  Bmc bmc_false(p_false, s);
  ProverResult r = bmc_false.check_until(2);
  EXPECT_EQ(r, ProverResult::FALSE);

  // add a reset
  Term reset_done = add_reset_seq(fts, rst, 1);
  // guard the property
  Term p_true_term = fts.make_term(Implies, reset_done, p_false_term);
  Property p(fts, p_true_term);
  // need to use a fresh solver to check the property again
  // pass the SolverEnum to use the same type of solver
  Bmc bmc(p, s->get_solver_enum());
  r = bmc.check_until(10);
  EXPECT_EQ(r,
            ProverResult::UNKNOWN);  // bmc can't prove, will only say unknown
}

TEST_P(ControlUnitTests, SimpleClock)
{
  // use a relational system -- easier to express positive clock edge without
  // syntactic functional constraints
  RelationalTransitionSystem rts(s);
  Term clk = rts.make_statevar("clk", boolsort);
  Term x = rts.make_statevar("x", bvsort8);
  // x' = (!clk & clk) ? x + 1 : x
  Term x_update =
      rts.make_term(Ite,
                    rts.make_term(And, rts.make_term(Not, clk), rts.next(clk)),
                    rts.make_term(BVAdd, x, rts.make_term(1, bvsort8)),
                    x);
  // should not be able to use assign next because the update is not functional
  // (contains next state var)
  EXPECT_THROW(rts.assign_next(x, x_update), PonoException);
  rts.constrain_trans(rts.make_term(Equal, x, x_update));
  rts.constrain_init(rts.make_term(Equal, x, rts.make_term(0, bvsort8)));
  Term state_counter = rts.make_statevar("state_counter", bvsort8);
  rts.assign_next(
      state_counter,
      rts.make_term(BVAdd, state_counter, rts.make_term(1, bvsort8)));
  rts.constrain_init(
      rts.make_term(Equal, state_counter, rts.make_term(0, bvsort8)));

  Term p_term = rts.make_term(
      Equal,
      x,
      rts.make_term(BVLshr, state_counter, rts.make_term(1, bvsort8)));
  Property p_false(rts, p_term);
  // use a new context so unroller doesn't clash on unrolled symbols
  Bmc bmc_false(p_false, s);
  ProverResult r = bmc_false.check_until(2);
  EXPECT_EQ(r, ProverResult::FALSE);

  // add a clock
  toggle_clock(rts, clk);
  // guard the property
  Property p(rts, p_term);
  // need to use a fresh solver to check the property again
  // pass the SolverEnum to use the same type of solver
  Bmc bmc(p, s->get_solver_enum());
  r = bmc.check_until(10);
  EXPECT_EQ(r,
            ProverResult::UNKNOWN);  // bmc can't prove, will only say unknown
}

INSTANTIATE_TEST_SUITE_P(ParameterizedControlUnitTests,
                         ControlUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
