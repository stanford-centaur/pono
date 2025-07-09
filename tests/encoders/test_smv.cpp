#include <string>
#include <tuple>

#include "core/rts.h"
#include "engines/kinduction.h"
#include "frontends/smv_encoder.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class SmvFileUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<
          tuple<SolverEnum, pair<const string, ProverResult>>>
{
};

TEST_P(SmvFileUnitTests, Encode)
{
  SmtSolver s = create_solver(get<0>(GetParam()));
  s->set_opt("incremental", "true");
  s->set_opt("produce-models", "true");
  RelationalTransitionSystem rts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  auto benchmark = get<1>(GetParam());
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/smv/";
  filename += benchmark.first;
  cout << "Reading file: " << filename << endl;
  SMVEncoder se(filename, rts);

  SafetyProperty prop(rts.solver(), se.propvec()[0]);
  KInduction kind(prop, rts, s);
  ProverResult res = kind.check_until(10);
  EXPECT_EQ(res, benchmark.second);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverSmvFileUnitTests,
    SmvFileUnitTests,
    testing::Combine(testing::ValuesIn(filter_solver_enums({ THEORY_INT })),
                     // from test_encoder_inputs.h
                     testing::ValuesIn(smv_inputs)));

}  // namespace pono_tests
