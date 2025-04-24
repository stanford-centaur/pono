#include <string>
#include <tuple>
#include <vector>

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
  SmtSolver s = create_solver(get<0>(GetParam()));
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/";
  filename += get<1>(GetParam());
  cout << "Reading file: " << filename << endl;
  BTOR2Encoder be(filename, fts);
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
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/mulo-test.btor2";
  BTOR2Encoder be(filename, fts);
  EXPECT_EQ(be.propvec().size(), 1);
  Property p(fts.solver(), be.propvec()[0]);
  KInduction kind(p, fts, s);
  ProverResult r = kind.check_until(2);
  EXPECT_EQ(r, ProverResult::TRUE);
}

TEST_P(Btor2UnitTests, InputConstraints)
{
  // test BTOR2 file with constraint containing input variables
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/mulo-test.btor2";
  BTOR2Encoder be(filename, fts);
  EXPECT_EQ(be.propvec().size(), 1);
  Property p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(6);
  ASSERT_NE(r, ProverResult::FALSE);
}

TEST_P(Btor2UnitTests, InputProp)
{
  // test BTOR2 file with bad containing input variables
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/input-in-bad.btor2";
  BTOR2Encoder be(filename, fts);
  EXPECT_EQ(be.propvec().size(), 1);
  Property p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(0);
  EXPECT_EQ(r, ProverResult::FALSE);
}

TEST_P(Btor2UnitTests, InvalidSmtlibSymbol)
{
  // test BTOR2 file with invalid SMT-LIB symbol
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/btor2/invalid-smtlib-symbol.btor2";
  BTOR2Encoder be(filename, fts);
  Property p(fts.solver(), be.propvec()[0]);
  Bmc bmc(p, fts, s);
  ProverResult r = bmc.check_until(0);
  ASSERT_NE(r, ProverResult::ERROR);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverBtor2FileUnitTests,
    Btor2FileUnitTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     // from test_encoder_inputs.h
                     testing::ValuesIn(btor2_inputs)));

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverBtor2UnitTests,
                         Btor2UnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests
