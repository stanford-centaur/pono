#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/rts.h"
#include "engines/bmc.h"
#include "engines/kinduction.h"
#include "modifiers/fault_injector.h"
#include "modifiers/single_bit_flip_fault.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class FaultUnitTests : public ::testing::Test,
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

TEST_P(FaultUnitTests, DefaultFaultInjection)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort8);
  fts.constrain_init(fts.make_term(Equal, x, fts.make_term(0, bvsort8)));
  fts.assign_next(
      x,
      fts.make_term(Ite,
                    fts.make_term(BVUlt, x, fts.make_term(10, bvsort8)),
                    fts.make_term(BVAdd, x, fts.make_term(1, bvsort8)),
                    fts.make_term(0, bvsort8)));

  Property prop(fts, fts.make_term(BVUle, x, fts.make_term(10, bvsort8)));
  KInduction kind(prop, s);
  ProverResult res = kind.check_until(2);
  EXPECT_EQ(res, ProverResult::TRUE);

  // TODO: test with fault injector
  //       need to fix state variables first
  // FaultInjector fi(fts);
  // FunctionalTransitionSystem & faulty_fts = fi.faulty_transition_system();
}

INSTANTIATE_TEST_SUITE_P(ParameterizedFaultUnitTests,
                         FaultUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
