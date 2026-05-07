#ifdef WITH_SLANG

#include <string>

#include "core/fts.h"
#include "engines/bmc.h"
#include "frontends/systemverilog_encoder.h"
#include "gtest/gtest.h"
#include "modifiers/mod_ts_prop.h"
#include "smt-switch/utils.h"
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

  // Encode `file` (which must contain a single assert property), run
  // BMC up to 5 cycles, and check that the result matches `expected`.
  void check_bmc(const string & file,
                 ProverResult expected = ProverResult::FALSE)
  {
    SmtSolver s = create_solver(GetParam());
    s->set_opt("incremental", "true");
    s->set_opt("produce-models", "true");
    FunctionalTransitionSystem fts(s);
    SystemVerilogEncoder enc(sv_path(file), fts);
    ASSERT_EQ(enc.propvec().size(), 1u);
    Term prop_term = enc.propvec()[0];

    // Promote any input vars referenced by the property to state vars
    // (with no update), matching the preprocessing done by pono.cpp.
    // This is required because SafetyProver only accepts predicates
    // over current-state variables.
    TransitionSystem ts = fts;
    if (!ts.only_curr(prop_term) && ts.no_next(prop_term)) {
      UnorderedTermSet ivs_in_prop;
      get_free_symbolic_consts(prop_term, ivs_in_prop);
      ts = promote_inputvars(ts, ivs_in_prop);
    }

    SafetyProperty prop(ts.solver(), prop_term);
    Bmc bmc(prop, ts, s);
    EXPECT_EQ(bmc.check_until(5), expected);
  }
};

// ---------------------------------------------------------------------------
// Block-level constructs
// ---------------------------------------------------------------------------

// Module ports + always_ff with non-blocking assignment + binary +.
TEST_P(SVUnitTests, EncodeCounter) { check_bmc("counter.sv"); }

// `initial` block sets register's initial state.
TEST_P(SVUnitTests, InitialBlock) { check_bmc("initial_block.sv"); }

// Sized integer literals (hex, binary, signed negative).
TEST_P(SVUnitTests, SizedLiterals) { check_bmc("sized_literals.sv"); }

// Multiple `assert property` statements yield multiple props.
TEST_P(SVUnitTests, MultipleAssertions)
{
  SmtSolver s = create_solver(GetParam());
  FunctionalTransitionSystem fts(s);
  SystemVerilogEncoder enc(sv_path("multi_assert.sv"), fts);
  EXPECT_EQ(enc.propvec().size(), 3u);
}

// Hierarchical SystemVerilog: a parent module instantiates a child
// that drives one of the parent's wires.  If sub-modules are encoded,
// the child's `assign` constrains the wire and BMC cannot falsify the
// assertion in any finite bound.  EXPECTED TO FAIL with the current
// encoder, which skips Instance members and leaves the wire as a free
// input -- BMC then picks values that violate the property.  Update
// the expectation once hierarchical encoding lands.
TEST_P(SVUnitTests, HierarchicalModules)
{
  check_bmc("hierarchical.sv", ProverResult::UNKNOWN);
}

// ---------------------------------------------------------------------------
// Statement kinds
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, IfElse) { check_bmc("if_else.sv"); }
TEST_P(SVUnitTests, CaseStatement) { check_bmc("case_stmt.sv"); }

// ---------------------------------------------------------------------------
// Expression kinds
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, Ternary) { check_bmc("ternary.sv"); }
TEST_P(SVUnitTests, BitSelect) { check_bmc("bit_select.sv"); }
TEST_P(SVUnitTests, RangeSelect) { check_bmc("range_select.sv"); }
TEST_P(SVUnitTests, Concat) { check_bmc("concat.sv"); }
TEST_P(SVUnitTests, Replication) { check_bmc("replication.sv"); }

// ---------------------------------------------------------------------------
// Operators
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, BinaryArith) { check_bmc("binary_arith.sv"); }
TEST_P(SVUnitTests, BinaryBitwise) { check_bmc("binary_bitwise.sv"); }
TEST_P(SVUnitTests, BinaryCompare) { check_bmc("binary_compare.sv"); }
TEST_P(SVUnitTests, BinaryLogical) { check_bmc("binary_logical.sv"); }
TEST_P(SVUnitTests, Shift) { check_bmc("shift.sv"); }
TEST_P(SVUnitTests, UnaryNot) { check_bmc("unary_not.sv"); }
TEST_P(SVUnitTests, UnaryMinus) { check_bmc("unary_minus.sv"); }
TEST_P(SVUnitTests, Reduction) { check_bmc("reduction.sv"); }

// ---------------------------------------------------------------------------
// Combinational logic
// ---------------------------------------------------------------------------

TEST_P(SVUnitTests, ContinuousAssign)
{
  check_bmc("continuous_assign.sv");
}

TEST_P(SVUnitTests, AlwaysComb) { check_bmc("always_comb.sv"); }

TEST_P(SVUnitTests, LegacyAlwaysStar)
{
  check_bmc("legacy_always_comb.sv");
}

INSTANTIATE_TEST_SUITE_P(ParameterizedSolverSVUnitTests,
                         SVUnitTests,
                         testing::ValuesIn(available_solver_enums()));

}  // namespace pono_tests

#endif
