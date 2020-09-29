#include <string>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "available_solvers.h"
#include "test_encoder_inputs.h"

#include "core/fts.h"
#include "engines/kinduction.h"
#include "frontends/btor2_encoder.h"

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
  Property p(fts, be.propvec()[0]);
  KInduction kind(p, s);
  ProverResult r = kind.check_until(2);
  EXPECT_EQ(r, ProverResult::TRUE);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverBtor2FileUnitTests,
    Btor2FileUnitTests,
    testing::Combine(testing::ValuesIn(available_no_logging_solver_enums()),
                     // from test_encoder_inputs.h
                     testing::ValuesIn(btor2_inputs)));

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverBtor2UnitTests,
    Btor2UnitTests,
    testing::ValuesIn(available_no_logging_solver_enums()));

}  // namespace pono_tests
