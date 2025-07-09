#include "engines/ceg_prophecy_arrays.h"
#include "engines/cegar_values.h"
#include "engines/ic3ia.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"

// need mathsat for ic3ia
#ifdef WITH_MSAT

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

TEST(CegValues, SimpleSafe)
{
  set_global_logger_verbosity(1);

  SmtSolver s = create_solver(MSAT);
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
      make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(prop, rts, s);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST(CegValues, SimpleUnsafe)
{
  set_global_logger_verbosity(1);

  SmtSolver s = create_solver(MSAT);
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
      make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(prop, rts, s);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::FALSE);
}

}  // namespace pono_tests

#endif
