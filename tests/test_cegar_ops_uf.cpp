#include <tuple>

#include "core/fts.h"
#include "engines/prover.h"
#include "gtest/gtest.h"
#include "options/options.h"
#include "smt/available_solvers.h"
#include "utils/make_provers.h"
#include "utils/ts_analysis.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

unordered_set<Engine> get_cegar_ops_uf_engines()
{
  return { BMC_SP, BMC, KIND };
}

unordered_set<Engine> get_cegar_ops_uf_bv_engines()
{
  return { BMC_SP, BMC, KIND, IC3SA_ENGINE };
}
unordered_set<Engine> get_cegar_ops_uf_interp_engines()
{
  return { DAR, IC3IA_ENGINE, INTERP, ISMC };
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

Term safe_prop(SmtSolver solver, const Term & x)
{
  auto sort = x->get_sort();
  auto lt_op = sort->get_sort_kind() == BV ? BVUlt : Lt;
  return solver->make_term(lt_op, x, solver->make_term(11, sort));
}

Term unsafe_prop(SmtSolver solver, const Term & x)
{
  auto sort = x->get_sort();
  auto lt_op = sort->get_sort_kind() == BV ? BVUlt : Lt;
  return solver->make_term(lt_op, x, solver->make_term(10, sort));
}

class CegOpsUfTestBase
{
 public:
  template <SortKind SK>
  void initialize_ts(SolverEnum se);

 protected:
  SmtSolver solver;
  Sort sort;
  Term x;
  FunctionalTransitionSystem fts;
};

template <>
void CegOpsUfTestBase::initialize_ts<BV>(SolverEnum se)
{
  solver = create_solver(se, se == BTOR);
  solver->set_opt("produce-unsat-assumptions", "true");
  sort = solver->make_sort(BV, 8);
  x = solver->make_symbol("x", sort);
  fts = counter_ts(solver, x);
}

template <>
void CegOpsUfTestBase::initialize_ts<INT>(SolverEnum se)
{
  // MathSAT does some rewriting that breaks invariant concretization.
  solver = create_solver(se, se == MSAT);
  solver->set_opt("produce-unsat-assumptions", "true");
  sort = solver->make_sort(INT);
  x = solver->make_symbol("x", sort);
  fts = counter_ts(solver, x);
}

template <SortKind SK>
class CegOpsUfTest
    : public testing::Test,
      public testing::WithParamInterface<tuple<Engine, bool, SolverEnum>>,
      public CegOpsUfTestBase
{
 protected:
  void SetUp() override
  {
    opts.engine_ = get<0>(GetParam());
    opts.ceg_bv_arith_as_free_symbol_ = get<1>(GetParam());
    opts.smt_solver_ = get<2>(GetParam());
    opts.check_invar_ = true;
    initialize_ts<SK>(opts.smt_solver_);
  }
  PonoOptions opts;
};

typedef CegOpsUfTest<BV> BVCegOpsUfTest;
typedef CegOpsUfTest<INT> IntCegOpsUfTest;

TEST_P(BVCegOpsUfTest, Safe)
{
  Term prop_term = safe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { BVAdd });
  ProverResult r = ceg_prover->check_until(20);
  if (opts.engine_ == Engine::BMC) {
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  } else {
    ASSERT_EQ(r, ProverResult::TRUE);
    if (opts.engine_ != Engine::BMC_SP && opts.engine_ != Engine::KIND) {
      Term invar = ceg_prover->invar();
      ASSERT_TRUE(check_invar(fts, prop_term, invar));
    }
  }
}

TEST_P(BVCegOpsUfTest, Unsafe)
{
  Term prop_term = unsafe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { BVAdd });
  ProverResult r = ceg_prover->check_until(10);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(ceg_prover->witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedCegOpsUfTests,
    BVCegOpsUfTest,
    testing::Combine(testing::ValuesIn(get_cegar_ops_uf_bv_engines()),
                     testing::ValuesIn({ false, true }),
                     testing::ValuesIn(filter_solver_enums({ THEORY_BV }))));

TEST_P(IntCegOpsUfTest, Safe)
{
  Term prop_term = safe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { Plus });
  ProverResult r = ceg_prover->check_until(15);
  if (opts.engine_ == Engine::BMC) {
    ASSERT_EQ(r, ProverResult::UNKNOWN);
  } else {
    ASSERT_EQ(r, ProverResult::TRUE);
    if (opts.engine_ != Engine::BMC_SP && opts.engine_ != Engine::KIND) {
      Term invar = ceg_prover->invar();
      ASSERT_TRUE(check_invar(fts, prop_term, invar));
    }
  }
}

TEST_P(IntCegOpsUfTest, Unsafe)
{
  Term prop_term = unsafe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { Plus });
  ProverResult r = ceg_prover->check_until(10);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(ceg_prover->witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedCegOpsUfTests,
    IntCegOpsUfTest,
    testing::Combine(testing::ValuesIn(get_cegar_ops_uf_engines()),
                     testing::ValuesIn({ false, true }),
                     testing::ValuesIn(filter_solver_enums({ THEORY_INT }))));

template <SortKind SK>
class CegOpsUfInterpTest : public testing::Test,
                           public testing::WithParamInterface<
                               tuple<Engine, bool, SolverEnum, SolverEnum>>,
                           public CegOpsUfTestBase
{
 protected:
  void SetUp() override
  {
    opts.engine_ = get<0>(GetParam());
    opts.ceg_bv_arith_as_free_symbol_ = get<1>(GetParam());
    opts.smt_solver_ = get<2>(GetParam());
    opts.smt_interpolator_ = get<3>(GetParam());
    opts.check_invar_ = true;
    initialize_ts<SK>(opts.smt_solver_);
  }
  PonoOptions opts;
};

typedef CegOpsUfInterpTest<BV> BVCegOpsUfInterpTest;
typedef CegOpsUfInterpTest<INT> IntCegOpsUfInterpTest;

TEST_P(BVCegOpsUfInterpTest, Safe)
{
  Term prop_term = safe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { BVAdd });
  ProverResult r = ceg_prover->check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = ceg_prover->invar();
  ASSERT_TRUE(check_invar(fts, prop_term, invar));
}

TEST_P(BVCegOpsUfInterpTest, Unsafe)
{
  if (opts.smt_interpolator_ == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5 get-interpolant fails for this test";
  }
  Term prop_term = unsafe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { BVAdd });

  ProverResult r = ceg_prover->check_until(10);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(ceg_prover->witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedCegOpsUfTests,
    BVCegOpsUfInterpTest,
    testing::Combine(
        testing::ValuesIn(get_cegar_ops_uf_interp_engines()),
        testing::ValuesIn({ false, true }),
        testing::ValuesIn(filter_solver_enums({ THEORY_BV })),
        testing::ValuesIn(filter_interpolator_enums({ THEORY_BV }))));

TEST_P(IntCegOpsUfInterpTest, Safe)
{
  if (opts.smt_interpolator_ == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5 get-interpolant does not terminate for this test";
  }
  Term prop_term = safe_prop(solver, x);
  SafetyProperty prop(solver, prop_term);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { Plus });
  ProverResult r = ceg_prover->check_until(15);
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = ceg_prover->invar();
  ASSERT_TRUE(check_invar(fts, prop_term, invar));
}

TEST_P(IntCegOpsUfInterpTest, Unsafe)
{
  if (opts.smt_interpolator_ == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5 get-interpolant does not terminate for this test";
  }
  Term prop_term = unsafe_prop(solver, x);
  shared_ptr<SafetyProver> ceg_prover = make_cegar_bv_arith_prover(
      opts.engine_, { solver, prop_term }, fts, solver, opts, { Plus });
  ProverResult r = ceg_prover->check_until(10);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(ceg_prover->witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedCegOpsUfTests,
    IntCegOpsUfInterpTest,
    testing::Combine(
        testing::ValuesIn(get_cegar_ops_uf_interp_engines()),
        testing::ValuesIn({ false, true }),
        testing::ValuesIn(filter_solver_enums({ THEORY_INT })),
        testing::ValuesIn(filter_interpolator_enums({ THEORY_INT }))));

}  // namespace pono_tests
