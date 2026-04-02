#include <tuple>

#include "engines/ceg_prophecy_arrays.h"
#include "engines/cegar_values.h"
#include "engines/ic3ia.h"
#include "gtest/gtest.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt/available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class CegarValuesTest
    : public testing::TestWithParam<tuple<SolverEnum, SolverEnum>>
{
 protected:
  void SetUp() override
  {
    opts.smt_solver_ = get<0>(GetParam());
    opts.smt_interpolator_ = get<1>(GetParam());
    s = create_solver(opts.smt_solver_);
    s->set_opt("produce-unsat-assumptions", "true");
  }
  PonoOptions opts;
  SmtSolver s;
};

TEST_P(CegarValuesTest, SimpleSafe)
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

  // TODO create a make_ command for this
  shared_ptr<SafetyProver> ceg =
      make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(prop, rts, s, opts);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(CegarValuesTest, SimpleUnsafe)
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
                    // off by one in the update
                    // e.g. using <= instead of <
                    rts.make_term(Le, d, rts.make_term(200, intsort)),
                    rts.make_term(Store, a, i, d),
                    a));

  Term prop_term = rts.make_term(
      Lt, rts.make_term(Select, a, j), rts.make_term(200, intsort));
  SafetyProperty prop(s, prop_term);

  // TODO create a make_ command for this
  shared_ptr<SafetyProver> ceg =
      make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(prop, rts, s, opts);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::FALSE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedCegarValuesTest,
    CegarValuesTest,
    testing::Combine(
        testing::ValuesIn(filter_solver_enums({ THEORY_INT })),
        testing::ValuesIn(filter_interpolator_enums({ THEORY_INT }))));

}  // namespace pono_tests
