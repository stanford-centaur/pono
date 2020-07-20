#ifdef WITH_COREIR
#define STRHELPER(A) #A
#define STRFY(A) STRHELPER(A)

#include <string>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"

#include "available_solvers.h"
#include "core/rts.h"
#include "frontends/coreir_encoder.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class CoreIRUnitTests : public ::testing::Test,
                        public ::testing::WithParamInterface<tuple<SolverEnum, string>>
{};

TEST_P(CoreIRUnitTests, Encode)
{
  SmtSolver s = available_solvers().at(get<0>(GetParam()))(false);
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
  SmtSolver s = available_solvers().at(get<0>(GetParam()))(false);
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
                     testing::ValuesIn(vector<string>{ "counters.json",
                                                       "WrappedPE_nofloats.json",
                                                       "SimpleALU.json" })));

}  // namespace pono_tests
#endif
