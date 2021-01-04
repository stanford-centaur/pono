#include <utility>
#include <vector>

#include "core/fts.h"
#include "engines/ic3sa.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/ts_analysis.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

TransitionSystem simple_uvw_system(const SmtSolver & solver)
{
  // simple system from IC3SA paper (with larger bitwidth):
  // Model Checking of Verilog RTL Using IC3 with Syntax-Guided Abstraction
  //   -- Aman Goel, Karem Sakallah
  FunctionalTransitionSystem fts(solver);
  Sort bvsort8 = fts.make_sort(BV, 8);
  Term u = fts.make_statevar("u", bvsort8);
  Term v = fts.make_statevar("v", bvsort8);
  Term w = fts.make_statevar("w", bvsort8);

  Term one = fts.make_term(1, bvsort8);
  fts.constrain_init(fts.make_term(Equal, u, one));
  fts.constrain_init(fts.make_term(Equal, v, one));
  fts.constrain_init(fts.make_term(Equal, w, one));

  Term cond =
      fts.make_term(Or, fts.make_term(BVUlt, u, v), fts.make_term(BVUlt, v, w));
  Term ifb = fts.make_term(BVAdd, u, v);
  Term elseb = fts.make_term(BVAdd, v, one);

  fts.assign_next(u, fts.make_term(Ite, cond, ifb, elseb));
  fts.assign_next(v, fts.make_term(BVAdd, v, one));
  fts.assign_next(w, fts.make_term(BVAdd, v, one));

  return fts;
}

class IC3SAUnitTests : public ::testing::Test,
                       public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override { s = create_solver(GetParam()); }
  SmtSolver s;
};

TEST_P(IC3SAUnitTests, SimpleSystemSafe)
{
  TransitionSystem ts = simple_uvw_system(s);
  Term u = ts.named_terms().at("u");
  Term v = ts.named_terms().at("v");
  Term one = ts.make_term(1, ts.make_sort(BV, 8));

  Property p(ts, ts.make_term(Distinct, ts.make_term(BVAdd, u, v), one));

  IC3SA ic3sa(p, s);
  ProverResult r = ic3sa.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
  ASSERT_TRUE(check_invar(ts, p.prop(), ic3sa.invar()));
}

TEST_P(IC3SAUnitTests, SimpleSystemUnsafe)
{
  TransitionSystem ts = simple_uvw_system(s);
  Term u = ts.named_terms().at("u");
  Term v = ts.named_terms().at("v");

  Sort bvsort8 = ts.make_sort(BV, 8);

  Property p(
      ts,
      ts.make_term(
          Distinct, ts.make_term(BVAdd, u, v), ts.make_term(10, bvsort8)));

  IC3SA ic3sa(p, s);
  ProverResult r = ic3sa.check_until(10);
  ASSERT_EQ(r, ProverResult::FALSE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverIC3SAUnitTests,
                         IC3SAUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
