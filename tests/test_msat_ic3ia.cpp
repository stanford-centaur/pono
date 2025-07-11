#ifdef WITH_MSAT_IC3IA

#include <vector>

#include "core/rts.h"
#include "engines/msat_ic3ia.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/ts_analysis.h"

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

  SafetyProperty p(s, rts.make_term(Ge, c, rts.make_term(0, intsort)));

  MsatIC3IA msat_ic3ia(p, rts, s);
  ProverResult res = msat_ic3ia.prove();
  EXPECT_EQ(res, ProverResult::TRUE);

  Term invar = msat_ic3ia.invar();
  EXPECT_TRUE(check_invar(rts, p.prop(), invar));
}

TEST_P(MsatIC3IAUnitTests, IntCounterUnsafe)
{
  size_t counter_limit = 10;
  RelationalTransitionSystem rts(s);
  Term c = rts.make_statevar("c", intsort);

  rts.constrain_init(rts.make_term(Equal, c, rts.make_term(0, intsort)));
  rts.constrain_trans(rts.make_term(
      Implies,
      rts.make_term(Lt, c, rts.make_term(counter_limit, intsort)),
      rts.make_term(Equal,
                    rts.next(c),
                    rts.make_term(Plus, c, rts.make_term(1, intsort)))));

  SafetyProperty p(s, rts.make_term(Ge, c, rts.make_term(0, intsort)));

  MsatIC3IA msat_ic3ia(p, rts, s);
  ProverResult res = msat_ic3ia.prove();
  EXPECT_EQ(res, ProverResult::FALSE);

  vector<UnorderedTermMap> witness;
  ASSERT_TRUE(msat_ic3ia.witness(witness));
  EXPECT_EQ(witness.size(), counter_limit + 2);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverMsatIC3IAUnitTests,
                         MsatIC3IAUnitTests,
                         testing::ValuesIn({ MSAT }));
}  // namespace pono_tests

#endif
