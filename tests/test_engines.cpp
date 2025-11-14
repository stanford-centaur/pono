#include <tuple>
#include <vector>

#include "core/fts.h"
#include "core/prop.h"
#include "core/proverresult.h"
#include "core/rts.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/dual_approx_reach.h"
#include "engines/interp_seq_mc.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/ts_analysis.h"

using namespace pono;
using namespace smt;

namespace pono_tests {

enum TSEnum
{
  Functional,
  Relational
};

class EngineTest
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, TSEnum>>
{
 protected:
  void SetUp() override
  {
    std::tuple<SolverEnum, TSEnum> t = GetParam();
    se = std::get<0>(t);
    TSEnum ts_type = std::get<1>(t);
    if (ts_type == Functional) {
      ts = new FunctionalTransitionSystem();
    } else {
      ts = new RelationalTransitionSystem();
    }

    Sort bvsort8 = ts->make_sort(BV, 8);
    Term max_val = ts->make_term(7, bvsort8);

    // populate TS with basic resetting counter system
    counter_system(*ts, max_val);

    Term x = ts->named_terms().at("x");

    Term true_prop = ts->make_term(BVUle, x, ts->make_term(7, bvsort8));
    true_p = new SafetyProperty(ts->solver(), true_prop);

    Term false_prop = ts->make_term(BVUle, x, ts->make_term(6, bvsort8));
    false_p = new SafetyProperty(ts->solver(), false_prop);
  }
  SolverEnum se;
  TransitionSystem * ts;
  SafetyProperty * true_p;
  SafetyProperty * false_p;
};

TEST_P(EngineTest, BmcTrue)
{
  SmtSolver s = create_solver(se);
  Bmc b(*true_p, *ts, s);
  ProverResult r = b.check_until(20);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(EngineTest, BmcFalse)
{
  SmtSolver s = create_solver(se);
  Bmc b(*false_p, *ts, s);
  ProverResult r = b.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(b.witness(cex));
}

TEST_P(EngineTest, BmcSimplePathTrue)
{
  SmtSolver s = create_solver(se);
  BmcSimplePath bsp(*true_p, *ts, s);
  ProverResult r = bsp.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(EngineTest, BmcSimplePathFalse)
{
  SmtSolver s = create_solver(se);
  BmcSimplePath bsp(*false_p, *ts, s);
  ProverResult r = bsp.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(bsp.witness(cex));
}

TEST_P(EngineTest, KInductionTrue)
{
  SmtSolver s = create_solver(se);
  KInduction kind(*true_p, *ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(EngineTest, KInductionFalse)
{
  SmtSolver s = create_solver(se);
  KInduction kind(*false_p, *ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(kind.witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedEngineTest,
    EngineTest,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(std::vector<TSEnum>{ Functional,
                                                            Relational })));

class InterpEngineTest : public EngineTest
{
 protected:
  void SetUp() override
  {
    EngineTest::SetUp();
    // configure the interpolator in the options
    opts.smt_interpolator_ = se;
    // use Bitwuzla as the base solver
    s = create_solver(BZLA);
  }
  SmtSolver s;
  PonoOptions opts;
};

TEST_P(InterpEngineTest, InterpTrue)
{
  InterpolantMC itpmc(*true_p, *ts, s, opts);
  ProverResult r = itpmc.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);

  Term invar = itpmc.invar();
  ASSERT_TRUE(check_invar(*ts, true_p->prop(), invar));
}

TEST_P(InterpEngineTest, InterpFalse)
{
  if (opts.smt_interpolator_ == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5 interpolation fails for this case";
  }
  InterpolantMC itpmc(*false_p, *ts, s, opts);
  ProverResult r = itpmc.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(itpmc.witness(cex));
}

TEST_P(InterpEngineTest, IsmcTrue)
{
  InterpSeqMC ismc(*true_p, *ts, s, opts);
  ProverResult r = ismc.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = ismc.invar();
  ASSERT_TRUE(check_invar(*ts, true_p->prop(), invar));
}

TEST_P(InterpEngineTest, IsmcFalse)
{
  InterpSeqMC ismc(*false_p, *ts, s, opts);
  ProverResult r = ismc.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(ismc.witness(cex));
}

TEST_P(InterpEngineTest, DarTrue)
{
  DualApproxReach dar(*true_p, *ts, s, opts);
  ProverResult r = dar.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = dar.invar();
  ASSERT_TRUE(check_invar(*ts, true_p->prop(), invar));
}

TEST_P(InterpEngineTest, DarFalse)
{
  DualApproxReach dar(*false_p, *ts, s, opts);
  ProverResult r = dar.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(dar.witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedInterpEngineTest,
    InterpEngineTest,
    testing::Combine(testing::ValuesIn(available_interpolator_enums()),
                     testing::ValuesIn({ Functional, Relational })));

SafetyProperty * create_non_inductive_system(TransitionSystem * ts)
{
  SmtSolver s = ts->solver();
  Sort boolsort = s->make_sort(BOOL);
  Sort bvsort8 = s->make_sort(BV, 8);

  // Simple non-inductive system/property
  Term cfg = ts->make_statevar("cfg", boolsort);
  Term a = ts->make_inputvar("a", bvsort8);
  Term b = ts->make_inputvar("b", bvsort8);
  Term initstate = ts->make_statevar("initstate", boolsort);
  Term out = ts->make_statevar("out", bvsort8);
  ts->add_constraint(s->make_term(
      Equal,
      out,
      s->make_term(
          Ite, cfg, s->make_term(BVAdd, a, b), s->make_term(BVSub, a, b))));
  ts->constrain_init(initstate);
  ts->assign_next(initstate, s->make_term(false));
  ts->assign_next(
      cfg, s->make_term(Or, cfg, initstate));  // keeps the same configuration

  // since property can only use state variables, we need a witness
  Term prop = s->make_term(Implies,
                           s->make_term(Not, initstate),
                           s->make_term(Equal, out, s->make_term(BVAdd, a, b)));
  Term witness = ts->make_statevar("witness", boolsort);
  ts->constrain_init(witness);
  ts->assign_next(witness, prop);
  return new SafetyProperty(ts->solver(), witness);
}

class NonInductiveTest : public ::testing::Test,
                         public ::testing::WithParamInterface<TSEnum>
{
 protected:
  void SetUp() override
  {
    // use Bitwuzla as the solver
    s = create_solver(BZLA);
    if (GetParam() == Functional) {
      ts = new FunctionalTransitionSystem(s);
    } else {
      ts = new RelationalTransitionSystem(s);
    }
    true_p = create_non_inductive_system(ts);
  }
  SmtSolver s;
  TransitionSystem * ts;
  SafetyProperty * true_p;
  PonoOptions opts;
};

TEST_P(NonInductiveTest, BmcFail)
{
  Bmc b(*true_p, *ts, s);
  ProverResult r = b.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(NonInductiveTest, BmcSimplePathWin)
{
  BmcSimplePath bsp(*true_p, *ts, s);
  ProverResult r = bsp.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(NonInductiveTest, KInductionWin)
{
  KInduction kind(*true_p, *ts, s);
  ProverResult r = kind.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(ParametrizedNonInductiveTest,
                         NonInductiveTest,
                         testing::ValuesIn({ Functional, Relational }));

class InterpNonInductiveTest
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, TSEnum>>
{
 protected:
  void SetUp() override
  {
    // use Bitwuzla as the base solver
    s = create_solver(BZLA);

    std::tuple<SolverEnum, TSEnum> t = GetParam();
    opts.smt_interpolator_ = std::get<0>(t);
    TSEnum ts_type = std::get<1>(t);

    if (ts_type == Functional) {
      ts = new FunctionalTransitionSystem(s);
    } else {
      ts = new RelationalTransitionSystem(s);
    }

    true_p = create_non_inductive_system(ts);
  }
  SmtSolver s;
  TransitionSystem * ts;
  SafetyProperty * true_p;
  PonoOptions opts;
};

TEST_P(InterpNonInductiveTest, InterpWin)
{
  if (opts.smt_interpolator_ == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5 interpolation fails for this case";
  }
  InterpolantMC itpmc(*true_p, *ts, s, opts);
  ProverResult r = itpmc.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedInterpNonInductiveTest,
    InterpNonInductiveTest,
    testing::Combine(testing::ValuesIn(available_interpolator_enums()),
                     testing::ValuesIn({ Functional, Relational })));

std::vector<PonoOptions> get_interp_options()
{
  PonoOptions default_opts;
  PonoOptions interp_first_and_last_props;
  interp_first_and_last_props.interp_props_ = INTERP_FIRST_AND_LAST_PROPS;
  PonoOptions eager_unroll;
  eager_unroll.interp_eager_unroll_ = true;
  PonoOptions backward_interp;
  backward_interp.interp_backward_ = true;
  PonoOptions no_frontier_simp;
  no_frontier_simp.interp_frontier_set_simpl_ = false;
  return { default_opts,
           interp_first_and_last_props,
           eager_unroll,
           backward_interp,
           no_frontier_simp };
}

class InterpOptionsTest
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, PonoOptions>>
{
 protected:
  void SetUp() override
  {
    // use Bitwuzla as the base solver
    s = create_solver(BZLA);
    Sort bvsort8 = s->make_sort(BV, 8);
    max_val = s->make_term(10, bvsort8);
    fts = FunctionalTransitionSystem(s);
    counter_system(fts, max_val);

    std::tuple<SolverEnum, PonoOptions> t = GetParam();
    opts = std::get<1>(t);
    opts.smt_interpolator_ = std::get<0>(t);
  }
  SmtSolver s;
  Term max_val;                    // max value for counter system
  FunctionalTransitionSystem fts;  // counter system
  PonoOptions opts;
};

TEST_P(InterpOptionsTest, CounterSystemUnsafe)
{
  if (opts.smt_interpolator_ == CVC5_INTERPOLATOR) {
    GTEST_SKIP() << "cvc5 interpolation fails for this case";
  }
  Term x = fts.named_terms().at("x");

  Term prop_term = s->make_term(BVUlt, x, max_val);
  SafetyProperty p(s, prop_term);

  InterpolantMC interp_mc(p, fts, s, opts);
  ProverResult r = interp_mc.prove();
  ASSERT_EQ(r, ProverResult::FALSE);
  std::vector<UnorderedTermMap> cex;
  ASSERT_TRUE(interp_mc.witness(cex));
}

TEST_P(InterpOptionsTest, CounterSystemSafe)
{
  Term x = fts.named_terms().at("x");

  Term prop_term = s->make_term(BVUle, x, max_val);
  SafetyProperty p(s, prop_term);

  InterpolantMC interp_mc(p, fts, s, opts);
  ProverResult r = interp_mc.prove();
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = interp_mc.invar();
  ASSERT_TRUE(check_invar(fts, prop_term, invar));
}

INSTANTIATE_TEST_SUITE_P(
    ParametrizedInterpOptionsTest,
    InterpOptionsTest,
    testing::Combine(testing::ValuesIn(available_interpolator_enums()),
                     testing::ValuesIn(get_interp_options())));

}  // namespace pono_tests
