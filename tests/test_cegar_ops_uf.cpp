#include "engines/cegar_ops_uf.h"
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

TEST(CegOpsUf, SimpleSafe)
{
  SmtSolver s = create_solver(MSAT);
  RelationalTransitionSystem rts(s);

  Sort sort = s->make_sort(BV, 8);
  Term max_val = rts.make_term(10, sort);
  Term x = rts.make_statevar("x", sort);
  SortKind sk = sort->get_sort_kind();
  PrimOp plus_op = (sk == BV) ? BVAdd : Plus;
  PrimOp lt_op = (sk == BV) ? BVUlt : Lt;
  Term inc_term = rts.make_term(plus_op, x, rts.make_term(1, sort));
  Term zero = rts.make_term(0, sort);
  rts.assign_next(
                 x, rts.make_term(Ite, rts.make_term(lt_op, x, max_val), inc_term, zero));
  rts.set_init(rts.make_term(Equal, x, zero));


  Term prop_term = rts.make_term(BVUlt, x, rts.make_term(11, sort));
  Property prop(s, prop_term);

  // TODO create a make_ command for this
  shared_ptr<Prover> ceg = make_shared<CegarOpsUf<IC3IA>>(prop, rts, s);

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);

}

}  // namespace pono_tests

#endif
