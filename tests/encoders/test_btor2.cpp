#include <algorithm>
#include <string>
#include <tuple>

#include "core/fts.h"
#include "core/prop.h"
#include "core/proverresult.h"
#include "engines/bmc.h"
#include "engines/kinduction.h"
#include "engines/kliveness.h"
#include "frontends/btor2_encoder.h"
#include "gtest/gtest.h"
#include "modifiers/liveness_to_safety_translator.h"
#include "modifiers/static_coi.h"
#include "options/options.h"
#include "smt-switch/smt.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"
#include "utils/exceptions.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

/** @param name the input's name, without the extension they all share
 *  @return the path to read it from
 */
string input_path(const string & name)
{
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  return string(STRFY(PONO_SRC_DIR)) + "/tests/encoders/inputs/btor2/" + name
         + ".btor2";
}

/** @param solver_enum the solver to create
 *  @return the solver, set up the way encoding a Btor2 input needs
 */
SmtSolver make_solver(SolverEnum solver_enum)
{
  SmtSolver solver = create_solver(solver_enum);
  solver->set_opt("incremental", "true");
  return solver;
}

/** One Btor2 input encoded into a system of its own. Checking a liveness
 *  property adds variables to the system it is checked in, so a test checking
 *  more than one needs more than one of these.
 */
struct Encoded
{
  Encoded(SolverEnum solver_enum, const string & name)
      : solver(make_solver(solver_enum)), fts(solver), be(input_path(name), fts)
  {
  }

  SmtSolver solver;
  FunctionalTransitionSystem fts;
  BTOR2Encoder be;
};

/** The condition set a generalized-Buchi search is given for a property: its
 *  justice conditions plus the file's fairness constraints, matching how
 *  pono.cpp combines them.
 *  @param be the encoder holding the parsed properties
 *  @param index which justice property to take the conditions from
 *  @return the combined condition set
 */
TermVec all_conditions(const BTOR2Encoder & be, size_t index = 0)
{
  TermVec conditions = be.justicevec().at(index);
  conditions.insert(conditions.end(), be.fairvec().begin(), be.fairvec().end());
  return conditions;
}

/** Runs an engine on a safety property.
 *  @param prop the property to check
 *  @param ts the system it belongs to
 *  @param bound how far to unroll
 *  @return what the engine found
 */
template <class Prover>
ProverResult check(const SafetyProperty & prop,
                   TransitionSystem & ts,
                   int bound)
{
  Prover prover(prop, ts, ts.solver());
  return prover.check_until(bound);
}

/** Translates a liveness property to safety and checks it, which modifies the
 *  system it is given.
 *  @param input the encoded system, which the translation adds to
 *  @param conditions what has to hold infinitely often
 *  @param bound how far to unroll
 *  @return what the engine found
 */
template <class Prover>
ProverResult check_l2s(Encoded & input, const TermVec & conditions, int bound)
{
  Term prop_term =
      LivenessToSafetyTranslator{}.translate(input.fts, conditions);
  return check<Prover>(
      SafetyProperty(input.solver, prop_term), input.fts, bound);
}

/** Checks a liveness property with k-liveness, which counts how often a
 *  condition is observed rather than translating it to safety.
 *  @param input the encoded system
 *  @param conditions what has to hold infinitely often
 *  @param bound how far to unroll
 *  @param options settings for the inner prover
 *  @return what k-liveness found
 */
ProverResult check_klive(Encoded & input,
                         const TermVec & conditions,
                         int bound,
                         PonoOptions options = PonoOptions())
{
  LivenessProperty prop(input.solver, conditions);
  KLiveness kliveness(prop, input.fts, input.solver, options);
  return kliveness.check_until(bound);
}

class Btor2UnitTests : public ::testing::Test,
                       public ::testing::WithParamInterface<SolverEnum>
{
};

class Btor2FileUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<tuple<SolverEnum, string>>
{
};

TEST_P(Btor2FileUnitTests, Encode)
{
  Encoded input(get<0>(GetParam()), get<1>(GetParam()));
  // make sure that all inputs in bad and constraint have been promoted
  UnorderedTermSet free_vars;
  for (const auto & c : input.fts.constraints()) {
    get_free_symbolic_consts(c.first, free_vars);
  }
  for (const auto & p : input.be.propvec()) {
    get_free_symbolic_consts(p, free_vars);
  }
  int num_input =
      count_if(free_vars.begin(), free_vars.end(), [&input](const Term & v) {
        return input.fts.is_input_var(v);
      });
  EXPECT_EQ(num_input, 0);
}

TEST_P(Btor2UnitTests, OverflowEncoding)
{
  Encoded input(GetParam(), "mulo_test");
  EXPECT_EQ(input.be.propvec().size(), 1);
  SafetyProperty prop(input.solver, input.be.propvec()[0]);
  EXPECT_EQ(check<KInduction>(prop, input.fts, 2), ProverResult::TRUE);
}

TEST_P(Btor2UnitTests, InputConstraints)
{
  // test BTOR2 file with constraint containing input variables
  Encoded input(GetParam(), "mulo_test");
  EXPECT_EQ(input.be.propvec().size(), 1);
  SafetyProperty prop(input.solver, input.be.propvec()[0]);
  ASSERT_NE(check<Bmc>(prop, input.fts, 6), ProverResult::FALSE);
}

TEST_P(Btor2UnitTests, VariablesKeepTheNamesFromTheFile)
{
  // options that name a variable, such as --reset, look it up by the name the
  // file gave it rather than the one the encoder generates from the line
  Encoded input(GetParam(), "state2input");
  ASSERT_EQ(input.be.statesvec().size(), 2);
  ASSERT_EQ(input.be.inputsvec().size(), 1);
  EXPECT_EQ(input.fts.lookup("state2input"), input.be.statesvec()[0]);
  EXPECT_EQ(input.fts.lookup("actualstate"), input.be.statesvec()[1]);
  EXPECT_EQ(input.fts.lookup("in"), input.be.inputsvec()[0]);
}

TEST_P(Btor2UnitTests, InputProp)
{
  // test BTOR2 file with bad containing input variables
  Encoded input(GetParam(), "input_in_bad");
  EXPECT_EQ(input.be.propvec().size(), 1);
  SafetyProperty prop(input.solver, input.be.propvec()[0]);
  EXPECT_EQ(check<Bmc>(prop, input.fts, 0), ProverResult::FALSE);
}

TEST_P(Btor2UnitTests, InvalidSmtlibSymbol)
{
  // test BTOR2 file with invalid SMT-LIB symbol
  Encoded input(GetParam(), "invalid_smtlib_symbol");
  SafetyProperty prop(input.solver, input.be.propvec()[0]);
  ASSERT_NE(check<Bmc>(prop, input.fts, 0), ProverResult::ERROR);
}

TEST_P(Btor2UnitTests, InitStateWithBool)
{
  // test BTOR2 file with bool init value
  Encoded input(GetParam(), "bool_init");
  SafetyProperty prop(input.solver, input.be.propvec()[0]);
  ASSERT_NE(check<Bmc>(prop, input.fts, 0), ProverResult::ERROR);
}

// A justice property with no fairness constraint: the lasso at mode=1
// satisfies the condition infinitely often, so it is a counterexample.
TEST_P(Btor2UnitTests, JusticeOnly)
{
  Encoded input(GetParam(), "justice_only");
  EXPECT_TRUE(input.be.fairvec().empty());
  EXPECT_EQ(all_conditions(input.be), input.be.justicevec().at(0));
  EXPECT_EQ(check_l2s<Bmc>(input, all_conditions(input.be), 10),
            ProverResult::FALSE);
}

// The same property under the other translator, which counts observations of
// the single condition rather than translating it to a safety property.
TEST_P(Btor2UnitTests, JusticeOnlyWithKLiveness)
{
  Encoded input(GetParam(), "justice_only");
  EXPECT_EQ(check_klive(input, all_conditions(input.be), 10),
            ProverResult::FALSE);
}

// A justice condition that holds once and never again: no lasso satisfies it
// infinitely often, so the property holds and refuting it is impossible.
TEST_P(Btor2UnitTests, JusticeHolds)
{
  Encoded input(GetParam(), "justice_holds");
  EXPECT_EQ(check_l2s<KInduction>(input, all_conditions(input.be), 20),
            ProverResult::TRUE);
}

// The same property under k-liveness, which proves it by bounding how often
// the condition can be observed. Its inner prover has to be one that can
// prove, so bmc, the default, would only ever come back unknown.
TEST_P(Btor2UnitTests, JusticeHoldsWithKLiveness)
{
  Encoded input(GetParam(), "justice_holds");
  PonoOptions options;
  options.engine_ = KIND;
  EXPECT_EQ(check_klive(input, all_conditions(input.be), 20, options),
            ProverResult::TRUE);
}

// A single justice property carrying two conditions, which is the generalized
// Buchi case the translator is written for. The toggling bit satisfies both
// infinitely often, so the property is violated.
TEST_P(Btor2UnitTests, JusticeWithTwoConditions)
{
  Encoded input(GetParam(), "justice_two_conditions");
  ASSERT_EQ(input.be.justicevec().at(0).size(), 2);
  EXPECT_EQ(check_l2s<Bmc>(input, all_conditions(input.be), 10),
            ProverResult::FALSE);
}

// Two conditions that cannot both recur, so the property holds. Checking them
// one at a time would refute it, which is what makes this worth pinning: the
// conditions have to be required together.
TEST_P(Btor2UnitTests, JusticeWithConflictingConditions)
{
  Encoded input(GetParam(), "justice_conflicting_conditions");
  ASSERT_EQ(input.be.justicevec().at(0).size(), 2);
  EXPECT_EQ(check_l2s<KInduction>(input, all_conditions(input.be), 20),
            ProverResult::TRUE);
}

// k-liveness counts observations of one condition, so a justice property with
// two of them is out of reach even with no fairness constraint in the file.
TEST_P(Btor2UnitTests, KLivenessRejectsTwoJusticeConditions)
{
  Encoded input(GetParam(), "justice_two_conditions");
  EXPECT_TRUE(input.be.fairvec().empty());
  LivenessProperty prop(input.solver, all_conditions(input.be));
  EXPECT_THROW(KLiveness(prop, input.fts, input.solver, PonoOptions()),
               PonoException);
}

// Each justice line is its own property, so a file can hold several and they
// need not agree. The second one here holds while the first is violated.
TEST_P(Btor2UnitTests, SecondJusticeProperty)
{
  Encoded input(GetParam(), "justice_two_properties");
  ASSERT_EQ(input.be.justicevec().size(), 2);
  EXPECT_EQ(check_l2s<KInduction>(input, all_conditions(input.be, 1), 20),
            ProverResult::TRUE);
}

// The cone of influence keeps what the justice conditions depend on and drops
// the rest, which has to leave the verdict alone. A state with an initial
// value could only be dropped from the transition relation, so this one has
// none.
TEST_P(Btor2UnitTests, StaticConeOfInfluenceKeepsJusticeState)
{
  Encoded input(GetParam(), "justice_unrelated_state");
  // The encoder names states after their line number, but lists them in the
  // order they are declared.
  ASSERT_EQ(input.be.statesvec().size(), 2);
  const Term mode = input.be.statesvec().at(0);
  const Term spare = input.be.statesvec().at(1);

  StaticConeOfInfluence coi(input.fts, all_conditions(input.be));
  const UnorderedTermSet & statevars = input.fts.statevars();
  EXPECT_EQ(statevars.size(), 1);
  EXPECT_TRUE(statevars.find(mode) != statevars.end());
  EXPECT_TRUE(statevars.find(spare) == statevars.end());

  EXPECT_EQ(check_l2s<Bmc>(input, all_conditions(input.be), 10),
            ProverResult::FALSE);
}

// The justice condition alone is violated by the lasso at mode=1, but that
// lasso never satisfies the fairness constraint, so adding the constraint
// removes the only counterexample. This is what makes the union of the two
// condition sets more than a no-op.
TEST_P(Btor2UnitTests, FairnessRestrictsCounterexample)
{
  const int bound = 10;

  Encoded justice(GetParam(), "fair_frozen_mode");
  EXPECT_EQ(check_l2s<Bmc>(justice, justice.be.justicevec().at(0), bound),
            ProverResult::FALSE);

  Encoded fair(GetParam(), "fair_frozen_mode");
  EXPECT_NE(check_l2s<Bmc>(fair, all_conditions(fair.be), bound),
            ProverResult::FALSE);

  // And no deeper counterexample exists either.
  Encoded deeper(GetParam(), "fair_frozen_mode");
  EXPECT_EQ(check_l2s<KInduction>(deeper, all_conditions(deeper.be), bound),
            ProverResult::TRUE);
}

// Fairness constraints that cannot all recur exclude every trace, so the
// property holds for want of any fair counterexample.
TEST_P(Btor2UnitTests, ContradictoryFairnessProvesVacuously)
{
  Encoded input(GetParam(), "fair_contradictory");
  ASSERT_EQ(input.be.fairvec().size(), 2);
  EXPECT_EQ(check_l2s<KInduction>(input, all_conditions(input.be), 10),
            ProverResult::TRUE);
}

// The cone has to be taken over the fairness constraints as well, not just
// the justice conditions. Taking it over the justice condition alone drops
// the state the constraint watches, and the translation then refers to a
// variable the system no longer has.
TEST_P(Btor2UnitTests, StaticConeOfInfluenceKeepsFairnessState)
{
  Encoded input(GetParam(), "fair_independent_state");
  ASSERT_EQ(input.be.statesvec().size(), 3);
  const Term mode = input.be.statesvec().at(0);
  const Term other = input.be.statesvec().at(1);
  const Term spare = input.be.statesvec().at(2);

  StaticConeOfInfluence coi(input.fts, all_conditions(input.be));
  const UnorderedTermSet & statevars = input.fts.statevars();
  EXPECT_EQ(statevars.size(), 2);
  EXPECT_TRUE(statevars.find(mode) != statevars.end());
  EXPECT_TRUE(statevars.find(other) != statevars.end());
  EXPECT_TRUE(statevars.find(spare) == statevars.end());

  EXPECT_EQ(check_l2s<Bmc>(input, all_conditions(input.be), 10),
            ProverResult::FALSE);
}

// k-liveness counts observations of a single condition, so it cannot honor a
// fairness constraint alongside a justice condition. It has to say so rather
// than check the justice condition on its own.
TEST_P(Btor2UnitTests, KLivenessRejectsFairness)
{
  Encoded input(GetParam(), "fair_frozen_mode");
  ASSERT_EQ(input.be.justicevec().at(0).size(), 1);
  ASSERT_EQ(input.be.fairvec().size(), 1);
  LivenessProperty prop(input.solver, all_conditions(input.be));
  try {
    KLiveness kliveness(prop, input.fts, input.solver, PonoOptions());
    FAIL() << "expected KLiveness to reject the fairness constraint";
  }
  catch (const PonoException & e) {
    EXPECT_NE(string(e.what()).find("fairness"), string::npos)
        << "message should name fairness constraints, got: " << e.what();
  }
}

INSTANTIATE_TEST_SUITE_P(
    ,
    Btor2FileUnitTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     // from test_encoder_inputs.h
                     testing::ValuesIn(btor2_inputs)),
    [](const auto & info) {
      return to_string(get<0>(info.param)) + "_" + get<1>(info.param);
    });

INSTANTIATE_TEST_SUITE_P(,
                         Btor2UnitTests,
                         testing::ValuesIn(available_solver_enums()),
                         testing::PrintToStringParamName());

}  // namespace pono_tests
