#ifdef WITH_COREIR

#include <string>
#include <tuple>
#include <vector>

#include "core/rts.h"
#include "frontends/coreir_encoder.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class CoreIRUnitTests
    : public ::testing::Test,
      public ::testing::WithParamInterface<tuple<SolverEnum, string>>
{
};

TEST_P(CoreIRUnitTests, Encode)
{
  SmtSolver s = create_solver(get<0>(GetParam()));
  RelationalTransitionSystem rts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/coreir/";
  filename += get<1>(GetParam());
  cout << "Reading file: " << filename << endl;
  CoreIREncoder ce(filename, rts);
}

TEST_P(CoreIRUnitTests, EncodeForceAbstract)
{
  SmtSolver s = create_solver(get<0>(GetParam()));
  RelationalTransitionSystem rts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/coreir/";
  filename += get<1>(GetParam());
  cout << "Reading file: " << filename << endl;
  CoreIREncoder ce(filename, rts, true);
}

INSTANTIATE_TEST_SUITE_P(
    ParameterizedSolverCoreIRUnitTests,
    CoreIRUnitTests,
    testing::Combine(testing::ValuesIn(available_solver_enums()),
                     // from test_encoder_inputs.h
                     testing::ValuesIn(coreir_inputs)));

}  // namespace pono_tests
#endif
