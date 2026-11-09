#include <string>
#include <tuple>

#include "core/tts.h"
#include "engines/kinduction.h"
#include "frontends/timed_vmt_encoder.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"
#include "utils/logger.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {
class TimedVMTFileUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<
          tuple<SolverEnum, pair<const string, ProverResult>>>
{
};

TEST_P(TimedVMTFileUnitTests, Encode)
{
  SmtSolver s = create_solver(get<0>(GetParam()));
  s->set_opt("incremental", "true");
  s->set_opt("produce-models", "true");
  TimedTransitionSystem tts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  auto benchmark = get<1>(GetParam());
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/tvmt/";
  filename += benchmark.first;
  TimedVMTEncoder tve(filename, tts);
  EXPECT_TRUE(tve.propvec().size() > 0);
  SafetyProperty prop(tts.solver(), tve.propvec()[0]);
  KInduction kind(prop, tts, s);
  ProverResult res = kind.check_until(10);
  EXPECT_EQ(res, benchmark.second);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverTimedVMTFileUnitTests,
    TimedVMTFileUnitTests,
    // testing::Combine(testing::ValuesIn({CVC5}),
    testing::Combine(testing::ValuesIn(filter_solver_enums({ THEORY_REAL })),
                     // from test_encoder_inputs.h
                     testing::ValuesIn(tvmt_inputs)));
}  // namespace pono_tests
