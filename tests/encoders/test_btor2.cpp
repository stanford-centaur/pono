#include <algorithm>
#include <string>
#include <tuple>

#include "core/fts.h"
#include "engines/bmc.h"
#include "engines/kinduction.h"
#include "frontends/btor2_encoder.h"
#include "gtest/gtest.h"
#include "smt-switch/utils.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"

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

class Btor2UnitTests : public ::testing::Test,
                       public ::testing::WithParamInterface<SolverEnum>
{
 protected:
  void SetUp() override
  {
    s = create_solver(GetParam());
    s->set_opt("incremental", "true");
  }

  SmtSolver s;
};

class Btor2FileUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<tuple<SolverEnum, string>>
{
 protected:
  void SetUp() override
  {
    s = create_solver(get<0>(GetParam()));
    s->set_opt("incremental", "true");
  }

  SmtSolver s;
};

TEST_P(Btor2FileUnitTests, Encode)
{
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path(get<1>(GetParam())), fts);
  // make sure that all inputs in bad and constraint have been promoted
  UnorderedTermSet free_vars;
  for (const auto & c : fts.constraints()) {
    get_free_symbolic_consts(c.first, free_vars);
  }
  for (const auto & p : be.propvec()) {
    get_free_symbolic_consts(p, free_vars);
  }
  int num_input =
      count_if(free_vars.begin(), free_vars.end(), [&fts](const Term & v) {
        return fts.is_input_var(v);
      });
  EXPECT_EQ(num_input, 0);
}

TEST_P(Btor2UnitTests, OverflowEncoding)
{
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("mulo_test"), fts);
  EXPECT_EQ(be.propvec().size(), 1);
  SafetyProperty p(fts.solver(), be.propvec()[0]);
  KInduction kind(p, fts, s);
  ProverResult r = kind.check_until(2);
  EXPECT_EQ(r, ProverResult::TRUE);
}

TEST_P(Btor2UnitTests, InputConstraints)
{
  // test BTOR2 file with constraint containing input variables
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("mulo_test"), fts);
  EXPECT_EQ(be.propvec().size(), 1);
  SafetyProperty p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(6);
  ASSERT_NE(r, ProverResult::FALSE);
}

TEST_P(Btor2UnitTests, VariablesKeepTheNamesFromTheFile)
{
  // options that name a variable, such as --reset, look it up by the name the
  // file gave it rather than the one the encoder generates from the line
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("state2input"), fts);
  ASSERT_EQ(be.statesvec().size(), 2);
  ASSERT_EQ(be.inputsvec().size(), 1);
  EXPECT_EQ(fts.lookup("state2input"), be.statesvec()[0]);
  EXPECT_EQ(fts.lookup("actualstate"), be.statesvec()[1]);
  EXPECT_EQ(fts.lookup("in"), be.inputsvec()[0]);
}

TEST_P(Btor2UnitTests, InputProp)
{
  // test BTOR2 file with bad containing input variables
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("input_in_bad"), fts);
  EXPECT_EQ(be.propvec().size(), 1);
  SafetyProperty p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(0);
  EXPECT_EQ(r, ProverResult::FALSE);
}

TEST_P(Btor2UnitTests, InvalidSmtlibSymbol)
{
  // test BTOR2 file with invalid SMT-LIB symbol
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("invalid_smtlib_symbol"), fts);
  SafetyProperty p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(0);
  ASSERT_NE(r, ProverResult::ERROR);
}

TEST_P(Btor2UnitTests, InitStateWithBool)
{
  // test BTOR2 file with bool init value
  FunctionalTransitionSystem fts(s);
  BTOR2Encoder be(input_path("bool_init"), fts);
  SafetyProperty p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(0);
  ASSERT_NE(r, ProverResult::ERROR);
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
