#include "core/fts.h"
#include "gtest/gtest.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class TSReplaceTests : public ::testing::Test,
                       public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    bvsort = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort bvsort;
};

TEST_P(TSReplaceTests, ReplaceTerms)
{
  FunctionalTransitionSystem fts(s);

  Term x = fts.make_statevar("x", bvsort);
  Term v = fts.make_statevar("v", bvsort);

  Term mul = s->make_term(BVMul, v, x);
  fts.assign_next(x, mul);

  Term fresh = fts.make_statevar("fresh", bvsort);
  fts.replace_terms({ { mul, fresh } });

  EXPECT_EQ(fresh, fts.state_updates().at(x));

  UnorderedTermSet free_syms;
  get_free_symbolic_consts(fts.init(), free_syms);
  get_free_symbolic_consts(fts.trans(), free_syms);

  EXPECT_TRUE(free_syms.find(v) == free_syms.end());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTSReplaceTests,
                         TSReplaceTests,
                         testing::ValuesIn(available_solver_enums()));
}  // namespace pono_tests
