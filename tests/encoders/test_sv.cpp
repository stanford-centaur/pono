#ifdef WITH_SLANG

#include <string>

#include "core/fts.h"
#include "engines/bmc.h"
#include "frontends/sv_encoder.h"
#include "gtest/gtest.h"
#include "smt/available_solvers.h"
#include "test_encoder_inputs.h"

using namespace pono;
using namespace smt;
using namespace std;

namespace pono_tests {

class SVUnitTests : public ::testing::Test,
                    public ::testing::WithParamInterface<SolverEnum>
{
};

TEST_P(SVUnitTests, EncodeCounter)
{
  SmtSolver s = create_solver(GetParam());
  s->set_opt("incremental", "true");
  s->set_opt("produce-models", "true");
  FunctionalTransitionSystem fts(s);
  // PONO_SRC_DIR is a macro set using CMake PROJECT_SRC_DIR
  string filename = STRFY(PONO_SRC_DIR);
  filename += "/tests/encoders/inputs/systemverilog/counter.sv";
  SVEncoder se(filename, fts);

  ASSERT_EQ(se.propvec().size(), 1);

  // counter.sv asserts (count != 4) on a 5-bit register that increments
  // every clock unless rst is high. With an unconstrained initial value
  // and a free rst input, BMC must find a counterexample within 5 steps.
  SafetyProperty prop(fts.solver(), se.propvec()[0]);
  Bmc bmc(prop, fts, s);
  ProverResult r = bmc.check_until(5);
  EXPECT_EQ(r, ProverResult::FALSE);
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnitTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
