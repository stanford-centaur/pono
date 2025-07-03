#include "engines/cegar_ops_uf.h"
#include "engines/ic3ia.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"

// need mathsat for ic3ia
#ifdef WITH_MSAT

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

RelationalTransitionSystem counter_ts(SmtSolver s, const Term & x)
{
  RelationalTransitionSystem rts(s);

  Sort sort = x->get_sort();
  Term nx = s->make_symbol("x.next", sort);
  rts.add_statevar(x, nx);

  Term max_val = rts.make_term(10, sort);
  SortKind sk = sort->get_sort_kind();
  PrimOp plus_op = (sk == BV) ? BVAdd : Plus;
  PrimOp lt_op = (sk == BV) ? BVUlt : Lt;
  Term inc_term = rts.make_term(plus_op, x, rts.make_term(1, sort));
  Term zero = rts.make_term(0, sort);

  rts.assign_next(
      x, rts.make_term(Ite, rts.make_term(lt_op, x, max_val), inc_term, zero));
  rts.set_init(rts.make_term(Equal, x, zero));

  return rts;
}

TEST(CegOpsUf, BVSimpleSafe)
{
  SmtSolver s = create_solver(MSAT);

  Sort sort = s->make_sort(BV, 8);
  Term x = s->make_symbol("x", sort);
  RelationalTransitionSystem rts = counter_ts(s, x);
  Term prop_term = rts.make_term(BVUlt, x, rts.make_term(11, sort));
  SafetyProperty prop(s, prop_term);

  shared_ptr<CegarOpsUf<IC3IA>> ceg =
      make_shared<CegarOpsUf<IC3IA>>(prop, rts, s);
  ceg->set_ops_to_abstract({ BVAdd });

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST(CegOpsUf, BVSimpleUnsafe)
{
  SmtSolver s = create_solver(MSAT);

  Sort sort = s->make_sort(BV, 8);
  Term x = s->make_symbol("x", sort);

  RelationalTransitionSystem rts = counter_ts(s, x);
  Term prop_term = rts.make_term(BVUlt, x, s->make_term(10, sort));
  SafetyProperty prop(s, prop_term);

  shared_ptr<CegarOpsUf<IC3IA>> ceg =
      make_shared<CegarOpsUf<IC3IA>>(prop, rts, s);
  ceg->set_ops_to_abstract({ BVAdd });

  ProverResult r = ceg->check_until(11);
  ASSERT_EQ(r, ProverResult::FALSE);
}

TEST(CegOpsUf, IntSimpleSafe)
{
  SmtSolver s = create_solver(MSAT);

  Sort sort = s->make_sort(INT);
  Term x = s->make_symbol("x", sort);
  RelationalTransitionSystem rts = counter_ts(s, x);
  Term prop_term = rts.make_term(Lt, x, rts.make_term(11, sort));
  SafetyProperty prop(s, prop_term);

  shared_ptr<CegarOpsUf<IC3IA>> ceg =
      make_shared<CegarOpsUf<IC3IA>>(prop, rts, s);
  ceg->set_ops_to_abstract({ Plus });

  ProverResult r = ceg->check_until(5);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST(CegOpsUf, IntSimpleUnsafe)
{
  SmtSolver s = create_solver(MSAT);

  Sort sort = s->make_sort(INT);
  Term x = s->make_symbol("x", sort);

  RelationalTransitionSystem rts = counter_ts(s, x);
  Term prop_term = rts.make_term(Lt, x, s->make_term(10, sort));
  SafetyProperty prop(s, prop_term);

  shared_ptr<CegarOpsUf<IC3IA>> ceg =
      make_shared<CegarOpsUf<IC3IA>>(prop, rts, s);
  ceg->set_ops_to_abstract({ Plus });

  ProverResult r = ceg->check_until(11);
  ASSERT_EQ(r, ProverResult::FALSE);
}

}  // namespace pono_tests

#endif
