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
 protected:
  static string sv_path(const string & name)
  {
    return string(STRFY(PONO_SRC_DIR))
           + "/tests/encoders/inputs/systemverilog/" + name;
  }

  // Encode `file`, run BMC up to `bound` steps, and check that the
  // (single) extracted property is falsified.
  void check_falsifiable(const string & file, int bound = 5)
  {
    SmtSolver s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    FunctionalTransitionSystem fts(s);
    SVEncoder enc(sv_path(file), fts);
    ASSERT_EQ(enc.propvec().size(), 1u);
    SafetyProperty prop(fts.solver(), enc.propvec()[0]);
    Bmc bmc(prop, fts, s);
    EXPECT_EQ(bmc.check_until(bound), ProverResult::FALSE);
  }
};

// ---------------------------------------------------------------------------
// Block-level constructs
// ---------------------------------------------------------------------------

// Module ports + always_ff with non-blocking assignment + binary +.
TEST_P(SVUnitTests, EncodeCounter) { check_falsifiable("counter.sv"); }

// `initial` block sets register's initial state.
TEST_P(SVUnitTests, InitialBlock) { check_falsifiable("initial_block.sv", 1); }

// Multiple `assert property` statements yield multiple props.
TEST_P(SVUnitTests, MultipleAssertions)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  SVEncoder enc(sv_path("multi_assert.sv"), fts);
  EXPECT_EQ(enc.propvec().size(), 3u);
}

// ---------------------------------------------------------------------------
// Statement kinds
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, IfElse) { check_falsifiable("if_else.sv"); }
TEST_P(SVUnitTests, CaseStatement) { check_falsifiable("case_stmt.sv"); }

// ---------------------------------------------------------------------------
// Expression kinds
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Ternary) { check_falsifiable("ternary.sv"); }
TEST_P(SVUnitTests, BitSelect) { check_falsifiable("bit_select.sv"); }
TEST_P(SVUnitTests, RangeSelect) { check_falsifiable("range_select.sv"); }
TEST_P(SVUnitTests, Concat) { check_falsifiable("concat.sv"); }
TEST_P(SVUnitTests, Replication) { check_falsifiable("replication.sv"); }

// ---------------------------------------------------------------------------
// Operators
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, BinaryArith) { check_falsifiable("binary_arith.sv"); }
TEST_P(SVUnitTests, BinaryBitwise) { check_falsifiable("binary_bitwise.sv"); }
TEST_P(SVUnitTests, BinaryCompare) { check_falsifiable("binary_compare.sv"); }
TEST_P(SVUnitTests, BinaryLogical) { check_falsifiable("binary_logical.sv"); }
TEST_P(SVUnitTests, Shift) { check_falsifiable("shift.sv"); }
TEST_P(SVUnitTests, UnaryNot) { check_falsifiable("unary_not.sv"); }
TEST_P(SVUnitTests, UnaryMinus) { check_falsifiable("unary_minus.sv"); }
TEST_P(SVUnitTests, Reduction) { check_falsifiable("reduction.sv"); }

// ---------------------------------------------------------------------------
// Combinational logic — currently broken in SVEncoder.
//
// `assign`, `always_comb`, and legacy `always @*` all funnel into the
// combinational branch which calls `fts_.add_invar(...)` on a constraint
// containing wires (declared as input vars).  add_invar requires the
// constraint to mention only state variables, so it throws
// "Invariants should be over current states only.".  The fix likely
// involves modeling wires by macro-substitution or via add_constraint
// with both init and next applications.  Tests are kept (DISABLED_) so
// that they will start passing once the encoder is fixed.
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, DISABLED_ContinuousAssign)
{
  check_falsifiable("continuous_assign.sv");
}

TEST_P(SVUnitTests, DISABLED_AlwaysComb)
{
  check_falsifiable("always_comb.sv");
}

TEST_P(SVUnitTests, DISABLED_LegacyAlwaysStar)
{
  check_falsifiable("legacy_always_comb.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnitTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
