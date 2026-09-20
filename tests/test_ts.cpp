#include "core/fts.h"
#include "core/prop.h"
#include "core/rts.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "utils/exceptions.h"
#include "utils/str_util.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class TSUnitTests : public ::testing::Test,
                    public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    // Boolector renames the symbols it is given, so it needs the logging
    // wrapper to report a variable under the name it was made with.
    s = create_solver(GetParam(), GetParam() == BTOR);
    bvsort = s->make_sort(BV, 8);
  }
  SmtSolver s;
  Sort bvsort;
};

TEST_P(TSUnitTests, FTS_IsFunc)
{
  FunctionalTransitionSystem fts(s);
  EXPECT_TRUE(fts.is_functional());

  // state variables without state updates
  // will make the system non-deterministic
  Term x = fts.make_statevar("x", bvsort);
  EXPECT_FALSE(fts.is_deterministic());
  EXPECT_TRUE(fts.is_functional());
  EXPECT_EQ(fts.statevars_with_no_update(), UnorderedTermSet({ x }));

  fts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  EXPECT_TRUE(fts.is_functional());
  EXPECT_TRUE(fts.is_deterministic());
  EXPECT_TRUE(fts.is_right_total());
  EXPECT_EQ(fts.statevars_with_no_update().size(), 0);

  fts.add_constraint(fts.make_term(BVUge, x, s->make_term(2, bvsort)));
  // any kind of constrains makes the system non-deterministic
  // TODO need to improve names here
  EXPECT_FALSE(fts.is_deterministic());
  EXPECT_FALSE(fts.is_right_total());

  Term y = fts.make_statevar("y", bvsort);
  fts.assign_next(y, y);
  EXPECT_TRUE(fts.is_functional());
  // still can't be deterministic because of the constraint
  EXPECT_FALSE(fts.is_deterministic());
  EXPECT_EQ(fts.statevars_with_no_update().size(), 0);

  Term ynext = fts.make_statevar("y.next", bvsort);
  EXPECT_EQ(fts.statevars_with_no_update(), UnorderedTermSet({ ynext }));

  TransitionSystem ts_copy = fts;
  EXPECT_EQ(fts.is_functional(), ts_copy.is_functional());
  EXPECT_EQ(fts.is_deterministic(), ts_copy.is_deterministic());
  EXPECT_EQ(fts.statevars_with_no_update(), ts_copy.statevars_with_no_update());
  EXPECT_EQ(ts_copy, fts);
  ts_copy.set_init(ts_copy.make_term(Equal, x, ts_copy.make_term(1, bvsort)));
  EXPECT_NE(ts_copy, fts);
}

TEST_P(TSUnitTests, RTS_IsFunc)
{
  RelationalTransitionSystem rts(s);
  EXPECT_FALSE(rts.is_functional());

  // state variables without state updates
  // will make the system non-functional
  Term x = rts.make_statevar("x", bvsort);
  EXPECT_FALSE(rts.is_functional());

  rts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  // Relational transition system is still not functional
  EXPECT_FALSE(rts.is_functional());
  // cannot guarantee determinism if relational
  EXPECT_FALSE(rts.is_deterministic());
  EXPECT_TRUE(rts.is_right_total());

  TransitionSystem ts_copy = rts;
  EXPECT_EQ(rts.is_functional(), ts_copy.is_functional());
  EXPECT_EQ(rts.is_deterministic(), ts_copy.is_deterministic());
  EXPECT_EQ(rts.is_right_total(), ts_copy.is_right_total());
  EXPECT_EQ(ts_copy, rts);
  ts_copy.set_init(ts_copy.make_term(Equal, x, ts_copy.make_term(1, bvsort)));
  EXPECT_NE(ts_copy, rts);
}

TEST_P(TSUnitTests, RTS_IsRightTotal_HierarchicalName)
{
  // Coverage test: is_right_total() builds quantifier param names by
  // concatenating ".param" onto Term::to_string(), which needs desanitizing
  // for names that require SMT-LIB `|...|` quoting (e.g. hierarchical names
  // like SystemVerilog produces, such as "mod.sig[3]"), matching the same
  // fix applied elsewhere. Unlike promote_inputvar() below, the available
  // backends don't validate quantifier param symbols strictly enough for
  // this to reproduce a failure pre-fix, but it still exercises the fixed
  // code path and checks the result is correct.
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("mod.sig[3]", bvsort);
  rts.assign_next(x, s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  bool right_total;
  ASSERT_NO_THROW(right_total = rts.is_right_total());
  EXPECT_TRUE(right_total);
}

TEST_P(TSUnitTests, PromoteInputvar_HierarchicalName)
{
  // Regression test: promote_inputvar() builds the next-state symbol's name
  // by concatenating a suffix onto Term::to_string(), which needs
  // desanitizing for names that require SMT-LIB `|...|` quoting.
  FunctionalTransitionSystem fts(s);
  Term iv = fts.make_inputvar("mod.sig[3]", bvsort);
  ASSERT_NO_THROW(fts.promote_inputvar(iv));
  EXPECT_TRUE(fts.is_curr_var(iv));
}

TEST_P(TSUnitTests, FTS_Exceptions)
{
  FunctionalTransitionSystem fts(s);
  Term x = fts.make_statevar("x", bvsort);
  Term xp1_n = fts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  EXPECT_THROW(fts.assign_next(x, xp1_n), PonoException);
}

TEST_P(TSUnitTests, RTS_Exceptions)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term xp1_n = rts.next(s->make_term(BVAdd, x, s->make_term(1, bvsort)));
  EXPECT_THROW(rts.assign_next(x, xp1_n), PonoException);
  EXPECT_NO_THROW(rts.constrain_trans(s->make_term(Equal, rts.next(x), xp1_n)));
}

TEST_P(TSUnitTests, FTS_DefaultCopy)
{
  FunctionalTransitionSystem fts;
  EXPECT_NO_THROW(fts = FunctionalTransitionSystem(s));
  // make sure the terms are not null
  EXPECT_TRUE(fts.init());
  EXPECT_TRUE(fts.trans());
}

TEST_P(TSUnitTests, RTS_Copy)
{
  RelationalTransitionSystem rts(s);

  RelationalTransitionSystem rts2 = rts;
  TransitionSystem ts = rts;
}

TEST_P(TSUnitTests, Prop_Copy)
{
  RelationalTransitionSystem rts(s);
  SafetyProperty p(s, s->make_term(true));

  SafetyProperty p2 = p;
}

TEST_P(TSUnitTests, RTS_ConstrainTrans)
{
  RelationalTransitionSystem rts(s);
  Term x = rts.make_statevar("x", bvsort);
  Term next_x = rts.next(x);
  rts.constrain_trans(rts.make_term(BVUge, x, next_x));
  EXPECT_TRUE(rts.statevars_with_no_update().empty());
}

// Pono marks the variables it adds for itself so the witness printers can
// tell them from the design's, which requires the marker to be on the name
// and the role the caller asked for to still be readable in it.
TEST_P(TSUnitTests, GeneratedVarsCarryTheMarker)
{
  RelationalTransitionSystem rts(s);
  const Term saved = rts.make_generated_statevar("saved", bvsort);
  EXPECT_TRUE(is_generated_name(name_desanitize(saved->to_string())));
  EXPECT_NE(saved->to_string().find("saved"), string::npos);

  const Term x = rts.make_statevar("x", bvsort);
  EXPECT_FALSE(is_generated_name(name_desanitize(x->to_string())));

  // A shadow of another variable names the variable it was derived from.
  const Term loop = rts.make_generated_statevar("loop", x, bvsort);
  EXPECT_TRUE(is_generated_name(name_desanitize(loop->to_string())));
  EXPECT_NE(loop->to_string().find("x"), string::npos);
}

// Asking twice for the same role is not a collision, since the caller cares
// about getting a variable rather than about the name it ends up with.
TEST_P(TSUnitTests, GeneratedVarsMoveAsideForEachOther)
{
  RelationalTransitionSystem rts(s);
  const Term first = rts.make_generated_statevar("monitor", bvsort);
  const Term second = rts.make_generated_statevar("monitor", bvsort);
  EXPECT_NE(first, second);
  EXPECT_TRUE(is_generated_name(name_desanitize(second->to_string())));
}

// A design variable named the way pono names its own would be left out of
// every trace, so it is refused rather than silently dropped.
TEST_P(TSUnitTests, DesignVarsCannotLookGenerated)
{
  RelationalTransitionSystem rts(s);
  const string taken = string(generated_prefix) + "saved";
  EXPECT_THROW(rts.make_statevar(taken, bvsort), PonoException);
  EXPECT_THROW(rts.make_inputvar(taken, bvsort), PonoException);

  // Frontends name a variable they made under a name of their own, which is
  // the other way a design's name reaches the system.
  const Term x = rts.make_statevar("x", bvsort);
  EXPECT_THROW(rts.name_term(taken, x), PonoException);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverTSUnitTests,
                         TSUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
