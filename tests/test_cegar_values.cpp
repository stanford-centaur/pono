#include "engines/ceg_prophecy_arrays.h"
#include "engines/cegar_values.h"
#include "engines/ic3ia.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/logger.h"
#include "utils/ts_analysis.h"

// need mathsat for ic3ia
#ifdef WITH_MSAT

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

TEST(CegValues, SimpleSafe)
{
  set_global_logger_verbosity(1);

  SmtSolver s = create_solver_for(MSAT, IC3IA_ENGINE, false, true);
  RelationalTransitionSystem rts(s);
  Sort bvsort8 = rts.make_sort(BV, 8);
  Sort bvsort32 = rts.make_sort(BV, 32);
  Sort arrsort = rts.make_sort(ARRAY, bvsort8, bvsort32);
  Term i = rts.make_statevar("i", bvsort8);
  Term j = rts.make_statevar("j", bvsort8);
  Term d = rts.make_statevar("d", bvsort32);
  Term a = rts.make_statevar("a", arrsort);

  Term constarr0 = rts.make_term(rts.make_term(0, bvsort32), arrsort);
  rts.set_init(rts.make_term(Equal, a, constarr0));
  rts.assign_next(
      a,
      rts.make_term(Ite,
                    rts.make_term(BVUlt, d, rts.make_term(200, bvsort32)),
                    rts.make_term(Store, a, i, d),
                    a));

  Term prop_term = rts.make_term(
      BVUlt, rts.make_term(Select, a, j), rts.make_term(200, bvsort32));
  Property prop(s, prop_term);

  // TODO create a make_ command for this
  shared_ptr<Prover> ceg =
      make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(prop, rts, s);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);

  Term invar = ceg->invar();
  cout << "got invariant " << invar << endl;
  // can't check invariant because it has prophecy variables in it
  // TODO consider checking the universally quantified version with CVC4
}

TEST(CegValues, SimpleUnsafe)
{
  set_global_logger_verbosity(1);

  SmtSolver s = create_solver_for(MSAT, IC3IA_ENGINE, false, true);
  RelationalTransitionSystem rts(s);
  Sort bvsort8 = rts.make_sort(BV, 8);
  Sort bvsort32 = rts.make_sort(BV, 32);
  Sort arrsort = rts.make_sort(ARRAY, bvsort8, bvsort32);
  Term i = rts.make_statevar("i", bvsort8);
  Term j = rts.make_statevar("j", bvsort8);
  Term d = rts.make_statevar("d", bvsort32);
  Term a = rts.make_statevar("a", arrsort);

  Term constarr0 = rts.make_term(rts.make_term(0, bvsort32), arrsort);
  rts.set_init(rts.make_term(Equal, a, constarr0));
  rts.assign_next(
      a,
      rts.make_term(Ite,
                    // off by one in the update
                    // e.g. using <= instead of <
                    rts.make_term(BVUle, d, rts.make_term(200, bvsort32)),
                    rts.make_term(Store, a, i, d),
                    a));

  Term prop_term = rts.make_term(
      BVUlt, rts.make_term(Select, a, j), rts.make_term(200, bvsort32));
  Property prop(s, prop_term);

  // TODO create a make_ command for this
  shared_ptr<Prover> ceg =
      make_shared<CegarValues<CegProphecyArrays<IC3IA>>>(prop, rts, s);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::FALSE);
}

}  // namespace pono_tests

#endif
