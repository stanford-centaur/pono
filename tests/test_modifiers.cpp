#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "smt-switch/utils.h"

#include "core/fts.h"
#include "core/rts.h"
#include "modifiers/history_refiner.h"
#include "modifiers/prophecy_refiner.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

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

TEST_P(ModifierUnitTests, HistoryRefiner)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  fts.constrain_init(fts.make_term(Equal, x, fts.make_term(0, bvsort)));
  fts.assign_next(x, fts.make_term(BVAdd, x, fts.make_term(1, bvsort)));

  HistoryRefiner hr(fts);

  size_t num_state_vars_orig = fts.statevars().size();

  Term trans_1 = fts.trans();
  Term hist_x_10 = hr.get_hist(x, 10);
  Term trans_2 = fts.trans();

  // should have added history variables to the transition relation
  EXPECT_NE(trans_1, trans_2);
  EXPECT_EQ(num_state_vars_orig + 10, fts.statevars().size());

  // shouldn't need to modify system
  // this one already needed to be created
  Term hist_x_2 = hr.get_hist(x, 2);
  Term trans_3 = fts.trans();

  EXPECT_EQ(trans_2, trans_3);
  EXPECT_EQ(num_state_vars_orig + 10, fts.statevars().size());
}

TEST_P(ModifierUnitTests, ProphecyRefinerSimple)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  fts.constrain_init(fts.make_term(Equal, x, fts.make_term(0, bvsort)));
  Term ite_update =
      fts.make_term(Ite,
                    fts.make_term(BVUlt, x, fts.make_term(9, bvsort)),
                    fts.make_term(BVAdd, x, fts.make_term(1, bvsort)),
                    x);
  fts.assign_next(x, ite_update);
  Term prop = fts.make_term(BVUlt, x, fts.make_term(10, bvsort));
  size_t num_statevars_orig = fts.statevars().size();
  EXPECT_EQ(num_statevars_orig, 1);

  ProphecyRefiner pr(fts);
  std::pair<Term, Term> p = pr.get_proph(x, 2, prop);
  Term proph_var = p.first;
  Term new_prop = p.second;

  TermVec free_vars;
  get_free_symbolic_consts(new_prop, free_vars);
  UnorderedTermSet free_vars_set(free_vars.begin(), free_vars.end());

  // Expecting two new history variable and one new prophecy variable
  EXPECT_EQ(fts.statevars().size() - 3, num_statevars_orig);

  // x should still be in the property
  EXPECT_TRUE(free_vars_set.find(x) != free_vars_set.end());

  // but now the prophecy variable should be also
  EXPECT_TRUE(free_vars_set.find(proph_var) != free_vars_set.end());
}

INSTANTIATE_TEST_SUITE_P(ParameterizedModifierUnitTests,
                         ModifierUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
