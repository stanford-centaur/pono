#include "core/fts.h"
#include "engines/prover.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/make_provers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

unordered_set<Engine> get_cegar_ops_uf_engines()
{
  return {
#ifdef WITH_MSAT
    IC3IA_ENGINE, INTERP,
#endif
    IC3SA_ENGINE, BMC_SP, BMC, KIND
  };
}

FunctionalTransitionSystem counter_ts(SmtSolver s, const Term & x)
{
  FunctionalTransitionSystem fts(s);

  Sort sort = x->get_sort();
  Term nx = s->make_symbol("x.next", sort);
  fts.add_statevar(x, nx);

  Term max_val = fts.make_term(10, sort);
  SortKind sk = sort->get_sort_kind();
  PrimOp plus_op = (sk == BV) ? BVAdd : Plus;
  PrimOp lt_op = (sk == BV) ? BVUlt : Lt;
  Term inc_term = fts.make_term(plus_op, x, fts.make_term(1, sort));
  Term zero = fts.make_term(0, sort);

  fts.assign_next(
      x, fts.make_term(Ite, fts.make_term(lt_op, x, max_val), inc_term, zero));
  fts.set_init(fts.make_term(Equal, x, zero));

  return fts;
}

class CegOpsUfTests : public ::testing::Test,
                      public ::testing::WithParamInterface<tuple<Engine, bool>>
{
 protected:
  void SetUp() override
  {
    tuple<Engine, bool> test_param = GetParam();
    opts.engine_ = get<0>(test_param);
    opts.ceg_bv_arith_as_free_symbol_ = get<1>(test_param);
    // use cvc5 as the base solver as it supports both BV and Int
    opts.smt_solver_ = SolverEnum::CVC5;
    solver = create_solver(opts.smt_solver_);
    solver->set_opt("produce-unsat-assumptions", "true");
  }
  PonoOptions opts;
  SmtSolver solver;
};

TEST_P(CegOpsUfTests, BVSimpleSafe)
{
  Sort sort = solver->make_sort(BV, 8);
  Term x = solver->make_symbol("x", sort);
  FunctionalTransitionSystem fts = counter_ts(solver, x);
  Term prop_term = fts.make_term(BVUlt, x, fts.make_term(11, sort));
  SafetyProperty prop(solver, prop_term);

  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, prop, fts, solver, opts, { BVAdd });

  ProverResult r = ceg_prover->check_until(5);
  if (opts.engine_ == Engine::BMC || opts.engine_ == Engine::BMC_SP) {
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  } else {
    ASSERT_EQ(r, ProverResult::TRUE);
  }
}

TEST_P(CegOpsUfTests, BVSimpleUnsafe)
{
  Sort sort = solver->make_sort(BV, 8);
  Term x = solver->make_symbol("x", sort);

  FunctionalTransitionSystem fts = counter_ts(solver, x);
  Term prop_term = fts.make_term(BVUlt, x, solver->make_term(10, sort));
  SafetyProperty prop(solver, prop_term);

  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, prop, fts, solver, opts, { BVAdd });

  ProverResult r = ceg_prover->check_until(11);
  ASSERT_EQ(r, ProverResult::FALSE);
}

TEST_P(CegOpsUfTests, IntSimpleSafe)
{
  if (opts.engine_ == Engine::IC3SA_ENGINE) {
    // IC3SA does not support Int
    return;
  }
  Sort sort = solver->make_sort(INT);
  Term x = solver->make_symbol("x", sort);
  FunctionalTransitionSystem fts = counter_ts(solver, x);
  Term prop_term = fts.make_term(Lt, x, fts.make_term(11, sort));
  SafetyProperty prop(solver, prop_term);

  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, prop, fts, solver, opts, { Plus });

  ProverResult r = ceg_prover->check_until(5);
  if (opts.engine_ == Engine::BMC || opts.engine_ == Engine::BMC_SP) {
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  } else {
    ASSERT_EQ(r, ProverResult::TRUE);
  }
}

TEST_P(CegOpsUfTests, IntSimpleUnsafe)
{
  if (opts.engine_ == Engine::IC3SA_ENGINE) {
    // IC3SA does not support Int
    return;
  }
  Sort sort = solver->make_sort(INT);
  Term x = solver->make_symbol("x", sort);

  FunctionalTransitionSystem fts = counter_ts(solver, x);
  Term prop_term = fts.make_term(Lt, x, solver->make_term(10, sort));
  SafetyProperty prop(solver, prop_term);

  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, prop, fts, solver, opts, { Plus });

  ProverResult r = ceg_prover->check_until(11);
  ASSERT_EQ(r, ProverResult::FALSE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedCegOpsUfTests,
    CegOpsUfTests,
    testing::Combine(testing::ValuesIn(get_cegar_ops_uf_engines()),
                     testing::ValuesIn({ false, true })));

}  // namespace pono_tests
