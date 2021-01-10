#include <utility>
#include <vector>

#include "core/fts.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "tests/common_ts.h"
#include "utils/exceptions.h"
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
    true_p = new Property(ts->solver(), true_prop);

    Term false_prop = ts->make_term(BVUle, x, ts->make_term(6, bvsort8));
    false_p = new Property(ts->solver(), false_prop);
  }
  SolverEnum se;
  Sort bvsort8;
  Term max_val;
  TransitionSystem * ts;
  Property * true_p;
  Property * false_p;
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
    true_p = new Property(ts->solver(), witness);

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
  Property * true_p;
};

TEST_P(InterpWinTests, BmcFail)
{
  Bmc b(*true_p, *ts, s);
  ProverResult r = b.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(InterpWinTests, BmcSimplePathFail)
{
  BmcSimplePath bsp(*true_p, *ts, s);
  ProverResult r = bsp.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(InterpWinTests, KInductionFail)
{
  KInduction kind(*true_p, *ts, s);
  ProverResult r = kind.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
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
#endif

}  // namespace pono_tests
