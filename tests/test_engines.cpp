#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/fts.h"
#include "core/rts.h"
#include "core/unroller.h"
#include "engines/bmc.h"
#include "engines/bmc_simplepath.h"
#include "engines/interpolantmc.h"
#include "engines/kinduction.h"
#include "utils/exceptions.h"

#include "available_solvers.h"

using namespace cosa;
using namespace smt;
using namespace std;

namespace cosa_tests {

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
    s = available_solvers().at(std::get<0>(t))();
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");

    TSEnum ts_type = std::get<1>(t);
    if (ts_type == Functional) {
      ts = new FunctionalTransitionSystem(s);
    } else {
      ts = new RelationalTransitionSystem(s);
    }

    bvsort8 = s->make_sort(BV, 8);

    // "Hello, World"-style counter test
    Term cnt = ts->make_state("cnt", bvsort8);
    ts->set_init(s->make_term(Equal, cnt, s->make_term(0, bvsort8)));
    ts->assign_next(
        cnt,
        s->make_term(Ite,
                     s->make_term(BVUle, cnt, s->make_term(6, bvsort8)),
                     s->make_term(BVAdd, cnt, s->make_term(1, bvsort8)),
                     s->make_term(0, bvsort8)));
    Term true_prop = s->make_term(BVUle, cnt, s->make_term(7, bvsort8));
    true_p = new Property(*ts, true_prop);

    Term false_prop = s->make_term(BVUle, cnt, s->make_term(6, bvsort8));
    false_p = new Property(*ts, false_prop);
  }
  SmtSolver s;
  Sort bvsort8;
  TransitionSystem * ts;
  Property * true_p;
  Property * false_p;
};

TEST_P(EngineUnitTests, BmcTrue)
{
  Bmc b(*true_p, s);
  ProverResult r = b.check_until(20);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(EngineUnitTests, BmcFalse)
{
  Bmc b(*false_p, s);
  ProverResult r = b.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
}

TEST_P(EngineUnitTests, BmcSimplePathTrue)
{
  BmcSimplePath bsp(*true_p, s);
  ProverResult r = bsp.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(EngineUnitTests, BmcSimplePathFalse)
{
  BmcSimplePath bsp(*false_p, s);
  ProverResult r = bsp.check_until(20);
  ASSERT_EQ(r, ProverResult::FALSE);
}

TEST_P(EngineUnitTests, KInductionTrue)
{
  KInduction kind(*true_p, s);
  ProverResult r = kind.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(EngineUnitTests, KInductionFalse)
{
  KInduction kind(*false_p, s);
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
    itp = ::smt::MsatSolverFactory::create_interpolating_solver();
  }
  SmtSolver itp;
};

TEST_P(InterpUnitTest, InterpTrue)
{
  InterpolantMC itpmc(*true_p, s, itp);
  ProverResult r = itpmc.check_until(20);
  ASSERT_EQ(r, ProverResult::TRUE);
}

TEST_P(InterpUnitTest, InterpFalse)
{
  InterpolantMC itpmc(*false_p, s, itp);
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
    s = ::smt::MsatSolverFactory::create();
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    itp = ::smt::MsatSolverFactory::create_interpolating_solver();

    TSEnum ts_type = GetParam();
    if (ts_type == Functional) {
      ts = new FunctionalTransitionSystem(s);
    } else {
      ts = new RelationalTransitionSystem(s);
    }

    boolsort = s->make_sort(BOOL);
    bvsort8 = s->make_sort(BV, 8);

    // Simple non-inductive system/property
    Term cfg = ts->make_state("cfg", boolsort);
    Term a = ts->make_input("a", bvsort8);
    Term b = ts->make_input("b", bvsort8);
    Term initstate = ts->make_state("initstate", boolsort);
    Term out = ts->make_state("out", bvsort8);
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
    Term witness = ts->make_state("witness", boolsort);
    ts->constrain_init(witness);
    ts->assign_next(witness, prop);
    true_p = new Property(*ts, witness);

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
  Bmc b(*true_p, s);
  ProverResult r = b.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(InterpWinTests, BmcSimplePathFail)
{
  BmcSimplePath bsp(*true_p, s);
  ProverResult r = bsp.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(InterpWinTests, KInductionFail)
{
  KInduction kind(*true_p, s);
  ProverResult r = kind.check_until(10);
  ASSERT_EQ(r, ProverResult::UNKNOWN);
}

TEST_P(InterpWinTests, InterpWin)
{
  InterpolantMC itpmc(*true_p, s, itp);
  ProverResult r = itpmc.check_until(10);
  ASSERT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedInterpWinTests,
                         InterpWinTests,
                         testing::ValuesIn({ Functional, Relational }));
#endif

}  // namespace cosa_tests
