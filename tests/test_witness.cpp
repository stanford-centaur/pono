#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class WitnessUnitTests : public ::testing::Test,
                         public ::testing::WithParamInterface<SolverEnum>
{
};

TEST_P(WitnessUnitTests, SimpleDefaultSolver)
{
  // use default solver
  FunctionalTransitionSystem fts;
  Sort bvsort8 = fts.make_sort(BV, 8);
  counter_system(fts, fts.make_term(20, bvsort8));
  Term x = fts.named_terms().at("x");

  Term eight = fts.make_term(8, bvsort8);
  Term prop_term = fts.make_term(BVUlt, x, eight);
  Property prop(fts.solver(), prop_term);

  SmtSolver s = create_solver(GetParam());
  Bmc bmc(prop, fts, s);
  ProverResult r = bmc.check_until(9);
  ASSERT_EQ(r, FALSE);

  vector<UnorderedTermMap> witness;
  bool ok = bmc.witness(witness);
  ASSERT_TRUE(ok);
  ASSERT_EQ(witness.size(), 9);
  ASSERT_NE(witness[8].find(x), witness[8].end());
  ASSERT_EQ(witness[8][x], eight);
}

TEST_P(WitnessUnitTests, ArraysDefaultSolver)
{
  // use default solver
  FunctionalTransitionSystem fts;
  Sort bvsort4 = fts.make_sort(BV, 4);
  Sort bvsort8 = fts.make_sort(BV, 8);
  Sort arrsort = fts.make_sort(ARRAY, bvsort4, bvsort8);

  Term arr = fts.make_statevar("arr", arrsort);
  Term idx = fts.make_statevar("idx", bvsort4);
  Term x = fts.make_statevar("x", bvsort4);
  Term data = fts.make_statevar("data", bvsort8);
  Term ext_idx = fts.make_term(Op(Zero_Extend, 4), idx);

  Term constarr0 = fts.make_term(fts.make_term(0, bvsort8), arrsort);
  fts.constrain_init(fts.make_term(Equal, arr, constarr0));
  fts.constrain_init(fts.make_term(Equal, idx, fts.make_term(0, bvsort4)));

  fts.assign_next(idx, fts.make_term(BVAdd, idx, fts.make_term(2, bvsort4)));
  fts.assign_next(arr, fts.make_term(Store, arr, idx, data));
  fts.add_invar(fts.make_term(BVUle, data, ext_idx));

  Term ten = fts.make_term(10, bvsort8);
  Term prop_term = fts.make_term(Distinct, fts.make_term(Select, arr, x), ten);
  Property prop(fts.solver(), prop_term);

  SmtSolver s = create_solver(GetParam());
  Bmc bmc(prop, fts, s);
  ProverResult r = bmc.check_until(6);
  ASSERT_EQ(r, FALSE);

  vector<UnorderedTermMap> witness;
  bool ok = bmc.witness(witness);
  ASSERT_TRUE(ok);
  ASSERT_EQ(witness.size(), 7);
  ASSERT_NE(witness[6].find(arr), witness[6].end());
  ASSERT_NE(witness[6].find(x), witness[6].end());
  ASSERT_EQ(witness[6][x], fts.make_term(10, bvsort4));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedWitnessUnitTests,
    WitnessUnitTests,
    testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
