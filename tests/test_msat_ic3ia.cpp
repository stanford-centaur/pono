#ifdef WITH_MSAT_IC3IA

#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "core/rts.h"
#include "engines/msat_ic3ia.h"
#include "utils/ts_analysis.h"

#include "available_solvers.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class MsatIC3IAUnitTests : public ::testing::Test,
                           public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    s->set_opt("produce-unsat-cores", "true");
    boolsort = s->make_sort(BOOL);
    intsort = s->make_sort(INT);
  }
  SmtSolver s;
  Sort boolsort, intsort;
};

TEST_P(MsatIC3IAUnitTests, IntCounterSafe)
{
  RelationalTransitionSystem rts(s);
  Term c = rts.make_statevar("c", intsort);

  rts.constrain_init(rts.make_term(Equal, c, rts.make_term(0, intsort)));
  rts.assign_next(c, rts.make_term(Plus, c, rts.make_term(1, intsort)));

  Property p(rts, rts.make_term(Ge, c, rts.make_term(0, intsort)));

  MsatIC3IA msat_ic3ia(p, s);
  ProverResult res = msat_ic3ia.prove();
  EXPECT_EQ(res, ProverResult::TRUE);
}

TEST_P(MsatIC3IAUnitTests, IntCounterUnsafe)
{
  RelationalTransitionSystem rts(s);
  Term c = rts.make_statevar("c", intsort);

  rts.constrain_init(rts.make_term(Equal, c, rts.make_term(0, intsort)));
  rts.constrain_trans(rts.make_term(Implies,
                                    rts.make_term(Lt, c, rts.make_term(10, intsort)),
                                    rts.make_term(Equal, rts.next(c), rts.make_term(Plus, c, rts.make_term(1, intsort)))));

  Property p(rts, rts.make_term(Ge, c, rts.make_term(0, intsort)));

  MsatIC3IA msat_ic3ia(p, s);
  ProverResult res = msat_ic3ia.prove();
  EXPECT_EQ(res, ProverResult::FALSE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverMsatIC3IAUnitTests,
                         MsatIC3IAUnitTests,
                         testing::ValuesIn({ MSAT }));
}  // namespace pono_tests

#endif
