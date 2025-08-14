#include <vector>

#include "core/fts.h"
#include "core/rts.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/dual_approx_reach.h"
#include "engines/interp_seq_mc.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/ts_analysis.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

enum TSEnum
{
  Functional,
  Relational
};

class EngineUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<std::tuple<SolverEnum, TSEnum>>
{
 protected:
  void SetUp() override
  {
    std::tuple<SolverEnum, TSEnum> t = GetParam();
    se = get<0>(t);
    TSEnum ts_type = std::get<1>(t);
    if (ts_type == Functional) {
      ts = new FunctionalTransitionSystem();
    } else {
      ts = new RelationalTransitionSystem();
    }

    bvsort8 = ts->make_sort(BV, 8);
    max_val = ts->make_term(7, bvsort8);

    // populate TS with basic resetting counter system
    counter_system(*ts, max_val);

    Term x = ts->named_terms().at("x");

    Term true_prop = ts->make_term(BVUle, x, ts->make_term(7, bvsort8));
    true_p = new SafetyProperty(ts->solver(), true_prop);

    Term false_prop = ts->make_term(BVUle, x, ts->make_term(6, bvsort8));
    false_p = new SafetyProperty(ts->solver(), false_prop);
  }
  SolverEnum se;
  Sort bvsort8;
  Term max_val;
  TransitionSystem * ts;
  SafetyProperty * true_p;
  SafetyProperty * false_p;
};

TEST_P(EngineUnitTests, BmcTrue)
{
  SmtSolver s = create_solver(se);
  Bmc b(*true_p, *ts, s);
  ProverResult r = b.check_until(20);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(EngineUnitTests, BmcFalse)
{
  SmtSolver s = create_solver(se);
  Bmc b(*false_p, *ts, s);
  ProverResult r = b.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(b.witness(cex));
}

TEST_P(EngineUnitTests, BmcSimplePathTrue)
{
  SmtSolver s = create_solver(se);
  BmcSimplePath bsp(*true_p, *ts, s);
  ProverResult r = bsp.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(EngineUnitTests, BmcSimplePathFalse)
{
  SmtSolver s = create_solver(se);
  BmcSimplePath bsp(*false_p, *ts, s);
  ProverResult r = bsp.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(bsp.witness(cex));
}

TEST_P(EngineUnitTests, KInductionTrue)
{
  SmtSolver s = create_solver(se);
  KInduction kind(*true_p, *ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(EngineUnitTests, KInductionFalse)
{
  SmtSolver s = create_solver(se);
  KInduction kind(*false_p, *ts, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(kind.witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedEngineUnitTests,
    EngineUnitTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     testing::ValuesIn(vector<TSEnum>{ Functional,
                                                       Relational })));

#if WITH_MSAT

class InterpUnitTest : public EngineUnitTests
{
 protected:
  void SetUp() override
  {
    EngineUnitTests::SetUp();
    s = create_solver(MSAT);
  }
  SmtSolver s;
};

TEST_P(InterpUnitTest, InterpTrue)
{
  InterpolantMC itpmc(*true_p, *ts, s);
  ProverResult r = itpmc.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);

  Term invar = itpmc.invar();
  ASSERT_TRUE(check_invar(*ts, true_p->prop(), invar));
}

TEST_P(InterpUnitTest, InterpFalse)
{
  InterpolantMC itpmc(*false_p, *ts, s);
  ProverResult r = itpmc.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(itpmc.witness(cex));
}

TEST_P(InterpUnitTest, IsmcTrue)
{
  InterpSeqMC ismc(*true_p, *ts, s);
  ProverResult r = ismc.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = ismc.invar();
  ASSERT_TRUE(check_invar(*ts, true_p->prop(), invar));
}

TEST_P(InterpUnitTest, IsmcFalse)
{
  InterpSeqMC ismc(*false_p, *ts, s);
  ProverResult r = ismc.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(ismc.witness(cex));
}

TEST_P(InterpUnitTest, DarTrue)
{
  DualApproxReach dar(*true_p, *ts, s);
  ProverResult r = dar.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
  Term invar = dar.invar();
  ASSERT_TRUE(check_invar(*ts, true_p->prop(), invar));
}

TEST_P(InterpUnitTest, DarFalse)
{
  DualApproxReach dar(*false_p, *ts, s);
  ProverResult r = dar.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(dar.witness(cex));
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedInterpUnitTest,
    InterpUnitTest,
    testing::Combine(testing::ValuesIn({ SolverEnum::MSAT }),
                     testing::ValuesIn({ Functional, Relational })));

class InterpWinTests : public ::testing::Test,
                       public ::testing::WithParamInterface<TSEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(MSAT);

    TSEnum ts_type = GetParam();
    if (ts_type == Functional) {
      ts = new FunctionalTransitionSystem(s);
    } else {
      ts = new RelationalTransitionSystem(s);
    }

    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);

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
    Term prop =
        s->make_term(Implies,
                     s->make_term(Not, initstate),
                     s->make_term(Equal, out, s->make_term(BVAdd, a, b)));
    Term witness = ts->make_statevar("witness", boolsort);
    ts->constrain_init(witness);
    ts->assign_next(witness, prop);
    true_p = new SafetyProperty(ts->solver(), witness);

    // debugging
    std::cout << "INIT" << std::endl;
    std::cout << ts->init() << std::endl;
    std::cout << "TRANS" << std::endl;
    std::cout << ts->trans() << std::endl;
    std::cout << "PROP" << std::endl;
    std::cout << witness << std::endl;
  }
  SmtSolver s;
  SmtSolver itp;
  Sort boolsort, bvsort8;
  TransitionSystem * ts;
  SafetyProperty * true_p;
};

TEST_P(InterpWinTests, BmcFail)
{
  Bmc b(*true_p, *ts, s);
  ProverResult r = b.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(InterpWinTests, BmcSimplePathWin)
{
  BmcSimplePath bsp(*true_p, *ts, s);
  ProverResult r = bsp.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(InterpWinTests, KInductionWin)
{
  KInduction kind(*true_p, *ts, s);
  ProverResult r = kind.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(InterpWinTests, InterpWin)
{
  InterpolantMC itpmc(*true_p, *ts, s);
  ProverResult r = itpmc.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedInterpWinTests,
                         InterpWinTests,
                         testing::ValuesIn({ Functional, Relational }));

vector<PonoOptions> get_interp_options()
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

class InterpOptionsTests : public ::testing::Test,
                           public ::testing::WithParamInterface<PonoOptions>
{
 protected:
  void SetUp() override
  {
    s = create_solver_for(MSAT, INTERP, false);
    bvsort8 = s->make_sort(BV, 8);
    max_val = s->make_term(10, bvsort8);
    fts = FunctionalTransitionSystem(s);
    counter_system(fts, max_val);
    opts = GetParam();
  }
  SmtSolver s;
  Sort bvsort8;
  Term max_val;                    // max value for counter system
  FunctionalTransitionSystem fts;  // counter system
  PonoOptions opts;
};

TEST_P(InterpOptionsTests, CounterSystemUnsafe)
{
  Term x = fts.named_terms().at("x");

  Term prop_term = s->make_term(BVUlt, x, max_val);
  SafetyProperty p(s, prop_term);

  InterpolantMC interp_mc(p, fts, s, opts);
  ProverResult r = interp_mc.prove();
  ASSERT_EQ(r, ProverResult::FALSE);
  vector<UnorderedTermMap> cex;
  ASSERT_TRUE(interp_mc.witness(cex));
}

TEST_P(InterpOptionsTests, CounterSystemSafe)
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

INSTANTIATE_TEST_SUITE_P(ParameterizedInterpOptionsTests,
                         InterpOptionsTests,
                         testing::ValuesIn(get_interp_options()));

#endif

}  // namespace pono_tests
