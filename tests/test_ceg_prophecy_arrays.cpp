#include "core/prop.h"
#include "core/rts.h"
#include "gtest/gtest.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt/available_solvers.h"
#include "utils/make_provers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class CegProphecyArraysTest : public ::testing::TestWithParam<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s->set_opt("produce-unsat-assumptions", "true");
  }
  SmtSolver s;
};

TEST_P(CegProphecyArraysTest, Simple)
{
  RelationalTransitionSystem rts(s);
  Sort intsort = rts.make_sort(INT);
  Sort arrsort = rts.make_sort(ARRAY, intsort, intsort);
  Term i = rts.make_statevar("i", intsort);
  Term j = rts.make_statevar("j", intsort);
  Term d = rts.make_statevar("d", intsort);
  Term a = rts.make_statevar("a", arrsort);

  Term constarr0 = rts.make_term(rts.make_term(0, intsort), arrsort);
  rts.set_init(rts.make_term(Equal, a, constarr0));
  rts.assign_next(
      a,
      rts.make_term(Ite,
                    rts.make_term(Lt, d, rts.make_term(200, intsort)),
                    rts.make_term(Store, a, i, d),
                    a));

  Term prop_term = rts.make_term(
      Lt, rts.make_term(Select, a, j), rts.make_term(200, intsort));
  SafetyProperty prop(s, prop_term);
  std::shared_ptr<SafetyProver> cegp =
      make_ceg_proph_prover(INTERP, prop, rts, s);
  ProverResult r = cegp->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedCegProphecyArraysTest,
    CegProphecyArraysTest,
    testing::ValuesIn(filter_solver_enums({ THEORY_INT })));

}  // namespace pono_tests
