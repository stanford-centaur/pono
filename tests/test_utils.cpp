#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/unroller.h"
#include "utils/exceptions.h"
#include "utils/make_provers.h"

#include "available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class UtilsUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, Engine>>
{
};

TEST_P(UtilsUnitTests, MakeProver)
{
  // use default solver
  FunctionalTransitionSystem fts;
  Sort bvsort8 = fts.make_sort(BV, 8);
  Sort boolsort = fts.make_sort(BOOL);
  Term one = fts.make_term(1, bvsort8);
  Term eight = fts.make_term(8, bvsort8);
  Term x = fts.make_statevar("x", bvsort8);

  fts.set_init(fts.make_term(Equal, x, fts.make_term(0, bvsort8)));
  fts.assign_next(x, fts.make_term(BVAdd, x, one));

  Term prop_term = fts.make_term(BVUlt, x, eight);
  Property prop(fts, prop_term);

  SolverEnum se = get<0>(GetParam());
  Engine eng = get<1>(GetParam());

  if (eng == INTERP && se != MSAT && se != MSAT_LOGGING) {
    // skip interpolation unless the solver is MathSAT
    return;
  }

  std::shared_ptr<Prover> prover = make_prover(eng, prop, se);
  ProverResult r = prover->check_until(9);
  ASSERT_EQ(r, FALSE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedUtilsUnitTests,
    UtilsUnitTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(all_engines())));

}  // namespace pono_tests
