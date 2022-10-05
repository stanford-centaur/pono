#include <utility>
#include <vector>

#include "core/rts.h"
#include "core/unroller.h"
#include "engines/kinduction.h"
#include "engines/mbic3.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class UFUnitTests : public ::testing::Test,
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

TEST_P(UFUnitTests, InductiveProp)
{
  // TODO: update this when boolector supports substitution for terms
  // with UF without logging
  if (s->get_solver_enum() == BTOR)

  {
    std::cout << "Warning: not running test with btor because it "
              << "doesn't support substitution (used in unrolling) for "
              << " terms containing UFs without using logging." << std::endl;
    return;
  }

  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term f = s->make_symbol("f", funsort);
  rts.constrain_init(rts.make_term(Equal, x, rts.make_term(1, bvsort)));
  rts.assign_next(x, rts.make_term(BVAdd, x, rts.make_term(1, bvsort)));
  rts.constrain_init(rts.make_term(Apply, f, rts.make_term(0, bvsort)));
  // f(x-1) -> f(x)
  // f(0) holds which is like the base case
  // because x starts at 1
  rts.constrain_trans(rts.make_term(
      Implies,
      rts.make_term(
          Apply, f, rts.make_term(BVSub, x, rts.make_term(1, bvsort))),
      rts.make_term(Apply, f, x)));

  Term p = rts.make_term(
      Apply, f, rts.make_term(BVSub, x, rts.make_term(1, bvsort)));
  Property prop(rts.solver(), p);

  s->push();
  KInduction kind(prop, rts, s);
  ProverResult r = kind.check_until(5);
  EXPECT_EQ(r, ProverResult::TRUE);
  s->pop();

  // TODO: re-enable support for UF in ModelBasedIC3
  // s->push();
  // ModelBasedIC3 ic3(prop, s);
  // r = ic3.check_until(10);
  // EXPECT_EQ(r, ProverResult::TRUE);
  // s->pop();
}

TEST_P(UFUnitTests, FalseProp)
{
  // TODO: update this when boolector supports substitution for terms
  // with UF without logging
  if (s->get_solver_enum() == BTOR) {
    std::cout << "Warning: not running btor because it doesn't support "
              << "substitution (used in unrolling) for terms containing "
              << "UFs without using logging." << std::endl;
    return;
  }

  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term f = s->make_symbol("f", funsort);
  rts.constrain_init(rts.make_term(Equal, x, rts.make_term(1, bvsort)));
  rts.assign_next(x, rts.make_term(BVAdd, x, rts.make_term(1, bvsort)));
  // this time not including the base case condition (makes it false)
  // rts.constrain_init(rts.make_term(Apply, f, rts.make_term(0, bvsort)));

  // f(x-1) -> f(x)
  // but this time f(0) isn't forced to hold
  rts.constrain_trans(rts.make_term(
      Implies,
      rts.make_term(
          Apply, f, rts.make_term(BVSub, x, rts.make_term(1, bvsort))),
      rts.make_term(Apply, f, x)));

  // guart property with a precondition so it doesn't fail in the initial state
  Term p = rts.make_term(
      Implies,
      rts.make_term(BVUge, x, rts.make_term(10, bvsort)),
      rts.make_term(
          Apply, f, rts.make_term(BVSub, x, rts.make_term(1, bvsort))));

  Property prop(rts.solver(), p);

  ASSERT_FALSE(rts.is_functional());
  s->push();
  KInduction kind(prop, rts, s);
  ProverResult r = kind.check_until(10);
  EXPECT_EQ(r, ProverResult::FALSE);
  s->pop();

  // TODO: re-enable support for UF in ModelBasedIC3 once we handle it correctly
  // s->push();
  // ModelBasedIC3 ic3(prop, s);
  // r = ic3.check_until(10);
  // EXPECT_EQ(r, ProverResult::FALSE);
  // s->pop();
}

INSTANTIATE_TEST_SUITE_P(ParameterizedUFUnitTests,
                         UFUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
