#include <utility>
#include <vector>

#include "core/fts.h"
#include "engines/ic3bits.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/ts_analysis.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class IC3BitsUnitTests : public ::testing::Test,
                         public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver_for(GetParam(), IC3_BITS, false);
    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort boolsort, bvsort8;
};

TEST_P(IC3BitsUnitTests, CounterSystemUnsafe)
{
  FunctionalTransitionSystem fts(s);
  Term max_val = fts.make_term(10, bvsort8);
  counter_system(fts, max_val);
  Term x = fts.named_terms().at("x");

  // off-by-one in property -- unsafe
  Term prop_term = s->make_term(BVUlt, x, max_val);
  Property p(s, prop_term);

  IC3Bits ic3bits(p, fts, s);
  ProverResult r = ic3bits.check_until(12);
  ASSERT_EQ(r, FALSE);
}

TEST_P(IC3BitsUnitTests, CounterSystemSafe)
{
  FunctionalTransitionSystem fts(s);
  Term max_val = fts.make_term(10, bvsort8);
  counter_system(fts, max_val);
  Term x = fts.named_terms().at("x");

  // safe property
  Term prop_term = s->make_term(BVUle, x, max_val);
  Property p(s, prop_term);

  IC3Bits ic3bits(p, fts, s);
  ProverResult r = ic3bits.prove();
  ASSERT_EQ(r, TRUE);
  Term invar = ic3bits.invar();
  ASSERT_TRUE(check_invar(fts, prop_term, invar));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverIC3BitsUnitTests,
    IC3BitsUnitTests,
    testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
